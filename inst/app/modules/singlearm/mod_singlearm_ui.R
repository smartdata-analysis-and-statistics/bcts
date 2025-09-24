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
            choices = c("Flat (Beta(1,1))" = "flat",
                        "Custom Beta prior" = "beta"),
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
          plotOutput(ns("sa_power_plot"))
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
      a_base <- if (prior_type == "beta") input$abase_sa else 1
      b_base <- if (prior_type == "beta") input$bbase_sa else 1

      B <- input$B_sa
      n_draws <- input$ndraws_sa
      seed <- input$seed_sa

      # Bayesian power
      power_res <- bcts::singlearm_beta_power(
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
      type1_res <- bcts::singlearm_beta_type1(
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

  })
}
