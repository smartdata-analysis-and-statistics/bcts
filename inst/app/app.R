# app.R
# Minimal Shiny app for Bayesian NI (Beta–Binomial, conjugate sampling)

library(shiny)
library(ggplot2)

# If your functions live in a package, uncomment:
# library(bcts)

ui <- fluidPage(
  withMathJax(),   # enable LaTeX rendering for the whole page
  titlePanel("Bayesian Trial Simulation (Beta–Binomial, conjugate)"),
  sidebarLayout(
    sidebarPanel(

      # --- Assumptions (truth used for simulation / planning) ---
      wellPanel(
        h4("Assumptions"),
        sliderInput("pc", HTML("True response rate in controls (\\(\\theta_c\\))"),
                    min = 0, max = 100, value = 85, step = 1,
                    post  = "%"),

        sliderInput("pt", HTML("True response rate in treatment (\\(\\theta_t\\))"),
                    min = 0, max = 100, value = 85, step = 1,
                    post  = "%"),

        selectInput(
          "prior",
          "Prior for control arm",
          choices = c("Flat (no external evidence)" = "flat", "Power prior (with historical data)" = "power"),
          selected = "flat"
        ),
        conditionalPanel(
          condition = "input.prior == 'power'",
          sliderInput("a0", "Discount factor a₀ (0 = ignore history, 1 = full borrow)",
                      min = 0, max = 100, value = 50, step = 1,post  = "%"),
          numericInput("y0", "Historical responders (y₀)", value = 64, min = 0, step = 1),
          numericInput("n0", "Historical sample size (n₀)", value = 75, min = 1, step = 1),
          numericInput("abase", "Baseline Beta a_base", value = 1, min = 0.001, step = 0.1),
          numericInput("bbase", "Baseline Beta b_base", value = 1, min = 0.001, step = 0.1)
        )
      ),

      # --- Design (what you plan to randomize) ---
      wellPanel(
        h4("Design"),
        numericInput("nc", "Number of patients in control arm",  value = 29, min = 1, step = 1),
        numericInput("nt", "Number of patients in treatment arm", value = 29, min = 1, step = 1),
        sliderInput(
          "M",
          label = HTML("Decision Margin (\\( \\Delta \\))"),
          min   = -100,
          max   =  100,
          value = -20,
          step  = 1,
          post  = "%"
        ),
        helpText(
          "Set Δ < 0 to evaluate non-inferiority (treatment may be up to |Δ| worse than control).",
          "Set Δ ≥ 0 to evaluate superiority (treatment must be at least Δ better). ",
          "Assumes higher response rates are better (responder events)."
        )
      ),



      # --- Posterior Evaluation ---
      wellPanel(
        h4("Posterior evaluation"),
        numericInput("ndraws", "Number of posterior draws", value = 2000, min = 1000, step = 100),
        helpText(HTML("Larger values give more precise results but increase runtime.")),
        checkboxInput("use_counts", "Use observed counts instead of simulating", value = FALSE),
        conditionalPanel(
          condition = "input.use_counts",
          numericInput("yc", "Observed y_c", value = NA, min = 0, step = 1),
          numericInput("yt", "Observed y_t", value = NA, min = 0, step = 1)
        ),
        numericInput("seed", "Seed (optional)", value = 123, step = 1),
      ),

      hr(),
      h4("Decision threshold & operating characteristics"),

      sliderInput("gamma", "Posterior threshold γ", min = 0.80, max = 0.99, value = 0.90, step = 0.001),
      checkboxInput("do_type1", "Estimate Type-I error at γ (LFN)", value = FALSE),
      checkboxInput("do_power", "Estimate Power at γ", value = FALSE),
      numericInput("B", "B (simulations for Type-I and Power)", value = 1000, min = 100, step = 100),

      actionButton("go", "Run")
    ),
    mainPanel(
      h4("Posterior NI result (single trial)"),
      verbatimTextOutput("single_out"),
      plotOutput("delta_plot", height = "320px"),
      hr(),
      h4("Operating characteristics at γ"),
      tableOutput("oc_table"),
      helpText("Type-I is evaluated at the least-favourable null p_t = p_c + M.")
    )
  )
)

server <- function(input, output, session) {

  # Build prior args list (reactive)
  prior_args <- reactive({
    if (input$prior == "flat") return(list())
    list(a0 = input$a0/100, y_0 = input$y0, n_0 = input$n0, a_base = input$abase, b_base = input$bbase)
  })

  # Run a single trial (or use observed y's) and compute posterior NI prob
  single <- eventReactive(input$go, {
    prior <- input$prior
    set.seed(if (!is.na(input$seed)) input$seed else NULL)

    if (isTRUE(input$use_counts) && !is.na(input$yc) && !is.na(input$yt)) {
      # Use observed counts; compute posterior shapes, then Pr(NI)
      sh <- posterior_shapes_conj(
        y_c = as.integer(input$yc),
        n_c = as.integer(input$nc),
        y_t = as.integer(input$yt),
        n_t = as.integer(input$nt),
        prior = prior,
        prior_args = prior_args()
      )
      prNI <- rbeta_diff_prob(
        a_t = sh$a_t, b_t = sh$b_t,
        a_c = sh$a_c, b_c = sh$b_c,
        M = input$M/100, n_draws = input$ndraws
      )
      list(
        mode = "observed",
        y_c = as.integer(input$yc), y_t = as.integer(input$yt),
        shapes = sh,
        prNI = prNI
      )
    } else {
      # Simulate a trial at (p_c, p_t, n_c, n_t)
      out <- bayesNI_trial_betaBinom_conj(
        p_c = input$pc/100, p_t = input$pt/100,
        n_c = as.integer(input$nc), n_t = as.integer(input$nt),
        M = input$M/100,
        prior = prior, prior_args = prior_args(),
        n_draws = input$ndraws,
        seed = if (!is.na(input$seed)) input$seed else NULL
      )
      csum <- out$summary
      list(
        mode = "simulated",
        y_c = as.integer(csum["y_c"]), y_t = as.integer(csum["y_t"]),
        shapes = out$shapes,
        prNI = as.numeric(csum["post_prob_NI"])
      )
    }
  }, ignoreInit = TRUE)

  # Print single-trial summary
  output$single_out <- renderPrint({
    s <- single()
    if (is.null(s)) return(invisible(NULL))
    cat(sprintf("Mode: %s\n", s$mode))
    cat(sprintf("y_c = %d / n_c = %d\n", s$y_c, input$nc))
    cat(sprintf("y_t = %d / n_t = %d\n", s$y_t, input$nt))
    cat(sprintf("Posterior Pr(Δ > %0.3f) = %0.4f\n", input$M/100, s$prNI))
    cat(sprintf("Decision at γ = %0.3f: %s\n",
                input$gamma, ifelse(s$prNI >= input$gamma, "Declare NI", "Do NOT declare NI")))
  })

  # Plot posterior of delta = theta_t - theta_c (MC from Beta posteriors)
  output$delta_plot <- renderPlot({
    s <- single()
    if (is.null(s)) return(invisible(NULL))

    # Draw posterior samples of theta_t and theta_c
    nd <- max(1000, input$ndraws)
    th_t <- stats::rbeta(nd, s$shapes$a_t, s$shapes$b_t)
    th_c <- stats::rbeta(nd, s$shapes$a_c, s$shapes$b_c)
    delta <- th_t - th_c

    ggplot(data.frame(delta = delta), aes(x = delta)) +
      geom_histogram(bins = 60) +
      geom_vline(xintercept = input$M/100, linetype = "dashed") +
      labs(x = expression(Delta == theta[t] - theta[c]),
           y = "Posterior density (MC histogram)",
           title = "Posterior of risk difference Δ",
           subtitle = sprintf("Pr(Δ > %0.3f) = %0.4f  |  γ = %0.3f",
                              input$M/100, s$prNI, input$gamma)) +
      theme_minimal(base_size = 12)
  })

  # Operating characteristics (optional)
  oc <- eventReactive(input$go, {
    res <- list()
    if (isTRUE(input$do_type1)) {
      res$type1 <- bayesNI_type1_betaBinom_conj(
        B = input$B,
        p_c = input$pc/100, M = input$M/100,
        n_c = as.integer(input$nc), n_t = as.integer(input$nt),
        threshold = input$gamma,
        prior = input$prior, prior_args = prior_args(),
        n_draws = input$ndraws,
        show_progress = TRUE
      )
    }
    if (isTRUE(input$do_power)) {
      res$power <- bayesNI_power_betaBinom_conj(
        B = input$B,
        p_c = input$pc/100, p_t = input$pt/100,
        n_c = as.integer(input$nc), n_t = as.integer(input$nt),
        M = input$M/100,
        threshold = input$gamma,
        prior = input$prior, prior_args = prior_args(),
        n_draws = input$ndraws,
        show_progress = TRUE
      )
    }
    res
  }, ignoreInit = TRUE)

  output$oc_table <- renderTable({
    r <- oc()
    if (length(r) == 0) return(NULL)
    out <- data.frame(
      Metric = character(0),
      Estimate = numeric(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(r$type1)) {
      out <- rbind(out, data.frame(Metric = "Type-I error", Estimate = r$type1$type1))
    }
    if (!is.null(r$power)) {
      out <- rbind(out, data.frame(Metric = "Power", Estimate = r$power))
    }
    out
  }, digits = 4)
}

shinyApp(ui, server)
