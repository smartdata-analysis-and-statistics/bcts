
#' @title Single-Arm Trial UI Module
#' @description UI for the single-arm Bayesian binomial trial tab
#' @export
mod_singlearm_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    title = "Single-arm trial",
    fluidPage(
      h4("Simulate a single-arm trial with a Bayesian decision rule"),
      p("This section will allow evaluation of single-arm designs using Beta–Binomial conjugate models."),
      p("You can simulate a posterior for a single group, compare against a threshold, or incorporate external data via power priors."),
      br(),
      sidebarLayout(
        sidebarPanel(
          h4("Design settings"),

          sliderInput(ns("pt_sa"), HTML("True response rate (\\( \\theta_t \\))"),
                      min = 0, max = 100, value = 80, step = 1, post = "%"),

          numericInput(ns("nt_sa"), HTML("Sample size (\\( n_t \\))"),
                       value = 40, min = 1, step = 1),

          selectInput(
            ns("prior_sa"),
            "Prior distribution",
            choices = c(
              "Flat (Beta(1,1))" = "flat",
              "Jeffreys (Beta(0.5,0.5))" = "jeffreys",
              "Custom Beta prior" = "beta"
            ),
            selected = "flat"
          ),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'beta'", ns("prior_sa")),
            numericInput(ns("abase_sa"), "Prior a_base", value = 1, min = 0.01, step = 0.1),
            numericInput(ns("bbase_sa"), "Prior b_base", value = 1, min = 0.01, step = 0.1)
          ),

          sat_decision_criteria_ui(ns("crit_sa"))
        ),

        mainPanel(
          h4("Design Summary"),
          htmlOutput(ns("sa_narrative_text")),
          verbatimTextOutput(ns("sa_summary")),

          br(),
          tabsetPanel(
            tabPanel("Bayesian Posterior", plotOutput(ns("bayes_posterior_plot"))),
            tabPanel("Pr(θ > M) vs y", plotOutput(ns("tail_prob_curve"))),
            tabPanel("Frequentist Test", plotOutput(ns("freq_binom_plot")))
          )
        )
      )
    )
  )
}


mod_singlearm_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    `%||%` <- function(x, y) if (is.null(x)) y else x
    fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
    fmt_int <- function(x) formatC(as.integer(x), big.mark = ",", format = "d")

    # Compute simulation results
    sa_results <- reactive({
      req(input$pt_sa, input$nt_sa,
          input[["crit_sa-M_sa"]],
          input[["crit_sa-gamma_sa"]])

      pt <- input$pt_sa / 100
      nt <- input$nt_sa
      M <- input[["crit_sa-M_sa"]] / 100
      gamma <- input[["crit_sa-gamma_sa"]] / 100

      prior_type <- input$prior_sa
      a_base <- switch(
        prior_type,
        "flat"     = 1,
        "jeffreys" = 0.5,
        "beta"     = input$abase_sa,
        1
      )
      b_base <- switch(
        prior_type,
        "flat"     = 1,
        "jeffreys" = 0.5,
        "beta"     = input$bbase_sa,
        1
      )

      B <- input$B_sa
      n_draws <- input$ndraws_sa
      seed <- input$seed_sa

      # Bayesian power
      power_res <- bcts::sat_betabinom_power(
        B = B,
        p_t = pt,
        n_t = nt,
        M = M,
        threshold = gamma,
        prior = prior_type,
        a_base = a_base,
        b_base = b_base,
        method = "exact",
        show_progress = FALSE
      )

      # Bayesian type-I error
      type1_res <- bcts::sat_betabinom_type1(
        B = B,
        n_t = nt,
        M = M,
        threshold = gamma,
        prior = prior_type,
        a_base = a_base,
        b_base = b_base,
        n_draws = n_draws,
        method = "exact",
        show_progress = FALSE
      )

      # Frequentist power
      crit_val <- qbinom(gamma, size = nt, prob = M) + 1
      freq_power <- 1 - pbinom(crit_val - 1, size = nt, prob = pt)


      # Frequentist type-I error
      freq_type1 <- 1 - gamma

      list(
        power = power_res,
        type1 = type1_res,
        freq_power = freq_power,
        freq_type1 = freq_type1
      )
    })


    # --- Render narrative directly ---
    output$sa_narrative_text <- renderText({
      pt <- input$pt_sa / 100
      nt <- input$nt_sa
      M  <- input[["crit_sa-M_sa"]] / 100
      gamma <- input[["crit_sa-gamma_sa"]] / 100
      alpha <- 0.05
      decision_mode <- "gamma"
      prior_type <- input$prior_sa
      a_base <- input$abase_sa %||% 1
      b_base <- input$bbase_sa %||% 1
      B <- input$B_sa
      ndraws <- input$ndraws_sa

      design_txt <- if (decision_mode == "gamma") {
        sprintf(
          "This single-arm design determines treatment efficacy by testing whether the posterior probability that the response rate exceeds %.0f%% is at least %.0f%%.",
          100 * M, 100 * gamma
        )
      } else {
        sprintf(
          "This design checks if the treatment response rate is above %.1f%% while keeping the chance of a false positive below %.1f%%.",
          100 * M, 100 * alpha
        )
      }

      arm_txt <- sprintf(
        "The trial includes %s patients, with an assumed true response rate of %.0f%%.",
        nt, 100 * pt
      )

      prior_txt <- switch(
        prior_type,
        "flat" = "The Bayesian analysis uses a flat prior, treating all response rates between 0% and 100% as equally likely before observing any data.",
        "beta" = sprintf(
          "The Bayesian analysis uses a Beta prior with shape parameters a = %.2f and b = %.2f, reflecting prior beliefs about likely response rates.",
          a_base, b_base
        ),
        "unknown prior type"
      )

      method_txt <- "Once data are observed, Bayesian inference is performed using a Beta–Binomial model. Because the Beta distribution is conjugate to the Binomial likelihood, the posterior distribution also follows a Beta distribution. This allows direct computation of the probability that the true response rate exceeds the pre-specified margin."

      HTML(paste(
        "<p>", design_txt, "</p>",
        "<p>", arm_txt, "</p>",
        "<p>", prior_txt, "</p>",
        "<p>", method_txt, "</p>"
      ))
    })

    # Text summary output
    output$sa_summary <- renderPrint({
      res <- sa_results()
      power_res <- res$power
      type1_res <- res$type1
      freq_power <- res$freq_power
      freq_type1 <- res$freq_type1

      cat("BAYESIAN ANALYSIS\n")
      cat(sprintf("Estimated power: %.2f%%\n", 100 * power_res$estimate))
      cat(sprintf("Estimated Type-I error: %.2f%%\n", 100 * type1_res$estimate))

      cat("\n")

      cat("FREQUENTIST COMPARISON\n")
      cat(sprintf("Frequentist power: %.2f%%\n", 100 * freq_power))
      cat(sprintf("Frequentist Type-I error: %.2f%%\n", 100 * freq_type1))
    })

    output$freq_binom_plot <- renderPlot({
      pt <- input$pt_sa / 100
      nt <- input$nt_sa
      M  <- input[["crit_sa-M_sa"]] / 100
      gamma <- input[["crit_sa-gamma_sa"]] / 100

      crit_val <- qbinom(gamma, size = nt, prob = M) + 1

      x_vals <- 0:nt
      probs <- dbinom(x_vals, size = nt, prob = M)
      barplot(probs, names.arg = x_vals, main = "Binomial Sampling Distribution",
              xlab = "Successes", ylab = "Probability",
              col = ifelse(x_vals >= crit_val, "red", "grey"))
      abline(v = crit_val, col = "red", lty = 2)
    })

    output$bayes_posterior_plot <- renderPlot({
      pt <- input$pt_sa / 100
      nt <- input$nt_sa
      M  <- input[["crit_sa-M_sa"]] / 100

      prior_type <- input$prior_sa
      a_base <- switch(prior_type,
                       "flat"     = 1,
                       "jeffreys" = 0.5,
                       "beta"     = input$abase_sa,
                       1)
      b_base <- switch(prior_type,
                       "flat"     = 1,
                       "jeffreys" = 0.5,
                       "beta"     = input$bbase_sa,
                       1)

      s <- round(pt * nt)
      a_post <- a_base + s
      b_post <- b_base + nt - s

      # Compute tail probability
      p_tail <- 1 - pbeta(M, a_post, b_post)
      p_tail_label <- sprintf("Pr(θ > M) = %s", sprintf("%.0f%%", 100 * p_tail))

      x_vals <- seq(0, 1, length.out = 500)
      df <- data.frame(
        theta = x_vals,
        prior = dbeta(x_vals, a_base, b_base),
        posterior = dbeta(x_vals, a_post, b_post)
      )

      ggplot(df, aes(x = theta)) +
        # Fill the posterior area
        geom_area(aes(y = posterior), fill = "firebrick", alpha = 0.4, show.legend = FALSE) +
        # Outline for posterior
        geom_line(aes(y = posterior, color = "Posterior"), linewidth = 1) +
        # Dashed prior line
        geom_line(aes(y = prior, color = "Prior"), linetype = "dashed", linewidth = 1) +
        # Tail shading beyond margin M
        geom_area(data = subset(df, theta >= M),
                  aes(y = posterior), fill = "firebrick", alpha = 0.2, show.legend = FALSE) +
        # Vertical line for M
        geom_vline(xintercept = M, color = "red", linetype = "dashed") +
        annotate("text", x = M, y = Inf, label = "M", vjust = 1.5, hjust = -0.2, color = "red") +
        # Add Pr(θ > M) annotation
        #annotate("text", x = M + 0.1, y = max(df$posterior) * 0.9,
        #         label = p_tail_label, color = "firebrick", size = 5, fontface = "italic") +
        scale_color_manual(values = c("Prior" = "grey40", "Posterior" = "firebrick")) +
        labs(
          title = "Prior and Posterior Distributions (Beta–Binomial)",
          x = expression(theta),
          y = "Density",
          fill = NULL,
          color = NULL,
          caption = sprintf(
            "Posterior is based on %d observed responses in %d patients and a Beta(%.1f, %.1f) prior, resulting in %s",
            s, nt, a_base, b_base, p_tail_label
          )
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.caption = element_text(size = 10, face = "italic", hjust = 0),
          legend.position = "top"
        )
    })

    output$tail_prob_curve <- renderPlot({
      pt <- input$pt_sa / 100  # treatment response rate
      nt <- input$nt_sa        # sample size
      M  <- input[["crit_sa-M_sa"]] / 100
      gamma <- input[["crit_sa-gamma_sa"]] / 100

      prior_type <- input$prior_sa
      a_base <- switch(prior_type,
                       "flat"     = 1,
                       "jeffreys" = 0.5,
                       "beta"     = input$abase_sa,
                       1)
      b_base <- switch(prior_type,
                       "flat"     = 1,
                       "jeffreys" = 0.5,
                       "beta"     = input$bbase_sa,
                       1)

      y_vals <- 0:nt

      # Compute tail probabilities for each possible y
      pr_theta_gt_M <- 1 - pbeta(M, a_base + y_vals, b_base + (nt - y_vals))
      exceeds_gamma <- pr_theta_gt_M > gamma

      # Frequentist likelihood under true pt
      likelihood_y <- dbinom(y_vals, size = nt, prob = pt)

      # Frequentist likelihood under the null hypothesis
      likelihood_y_type1  <- dbinom(y_vals, size = nt, prob = M)

      # Bayesian power: sum of probabilities where posterior exceeds gamma
      bayes_power <- sum(likelihood_y[exceeds_gamma])

      # Bayesian type1: sum of probabilities where posterior exceeds gamma
      bayes_type1 <- sum(likelihood_y_type1[exceeds_gamma])

      df_power <- data.frame(
        y = y_vals,
        pr_theta_gt_M,
        exceeds_gamma,
        prob = likelihood_y,
        metric = "Power"
      )

      df_type1 <- data.frame(
        y = y_vals,
        pr_theta_gt_M,
        exceeds_gamma,
        prob = likelihood_y_type1,
        metric = "Type-I Error"
      )

      df_plot <- rbind(df_power, df_type1)

      ggplot(df_plot, aes(x = y, y = pr_theta_gt_M)) +
        geom_col(aes(fill = exceeds_gamma, alpha = prob), width = 0.8) +
        geom_hline(yintercept = gamma, linetype = "dashed", color = "red", linewidth = 1) +
        scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey80")) +
        scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
        facet_wrap(~ metric, ncol = 1)  +
        labs(
          title = expression("Posterior probability " * Pr(theta > M) * " across outcomes"),
          subtitle = sprintf(
            "Power assumes true response rate = %.0f%%; Type-I Error assumes θ = M = %.0f%%",
            100 * pt, 100 * M
          ),
          caption = paste(
            "Posterior probabilities are computed across all possible outcomes y = 0, ..., n,",
            "using a Beta prior and Binomial likelihood.",
            "Bar opacity reflects the likelihood of observing y given the assumed response rate.",
            "Power is based on Pr(theta > M) under the true response rate (pt),",
            "while Type-I error assumes theta equals the decision threshold M.",
            sep = "\n"
          ),
          x = "Number of observed responses (y)",
          y = expression("Pr(" * theta * " > M)"),
          fill = NULL
        ) +
        theme_minimal(base_size = 13) +
        theme(
          legend.position = "none",
          plot.caption = element_text(hjust = 0),              # 0 = left, 0.5 = center, 1 = right
          plot.caption.position = "plot"                       # ensures caption aligns as part of the plot area
              )
    })

  })
}
