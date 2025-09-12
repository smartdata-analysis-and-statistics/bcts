# app.R
# Minimal Shiny app for Bayesian NI (Beta–Binomial, conjugate sampling)

library(shiny)
library(ggplot2)


# ---- small helper: get Pr(NI) across B trials for plotting ------------------
prNI_draws_conj <- function(B, p_c, p_t, n_c, n_t, M,
                            prior = c("flat","power"), prior_args = list(),
                            n_draws = 2000, seed = NULL) {
  prior <- match.arg(prior)
  if (!is.null(seed)) set.seed(seed)
  pr <- numeric(B)
  for (b in seq_len(B)) {
    res <- bayesNI_trial_betaBinom_conj(
      p_c = p_c, p_t = p_t, n_c = n_c, n_t = n_t, M = M,
      prior = prior, prior_args = prior_args,
      n_draws = n_draws
    )
    pr[b] <- unname(res$summary["post_prob_NI"])
  }
  pr
}


# If your functions live in a package, uncomment:
# library(bcts)

ui <- fluidPage(
  withMathJax(),   # enable LaTeX rendering for the whole page
  titlePanel("Bayesian Trial Simulation (Beta–Binomial, conjugate)"),
  sidebarLayout(
    sidebarPanel(

      # --- Assumptions (truth used for simulation / planning) ---
      wellPanel(
        h4("Treatment Arm Assumptions"),
        sliderInput("pt", HTML("True response rate  (\\(\\theta_t\\))"),
                    min = 0, max = 100, value = 85, step = 1,
                    post  = "%"),
        numericInput("nt", HTML("Number of randomized patients  (\\(n_t \\))"), value = 29, min = 1, step = 1),
      ),

      # --- Assumptions (truth used for simulation / planning) ---
      wellPanel(
        h4("Control Arm Assumptions"),
        sliderInput("pc", HTML("True response rate (\\(\\theta_c\\))"),
                    min = 0, max = 100, value = 85, step = 1,
                    post  = "%"),
        numericInput("nc", HTML("Number of randomized patients (\\(n_c \\))"),  value = 29, min = 1, step = 1),

        selectInput(
          "prior",
          "Prior distribution",
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


      wellPanel(
        h4("Decision criteria"),
        sliderInput(
          "M",
          label = HTML("Decision Margin (\\( \\Delta \\))"),
          min   = -100, max = 100, value = -20, step = 1, post = "%"
        ),
        helpText(
          "Δ < 0: non-inferiority (treatment may be up to |Δ| worse).",
          "Δ ≥ 0: superiority (treatment must be at least Δ better).",
          "Assumes higher response rates are better (responder events)."
        ),

        # Choice: set gamma directly OR set alpha and calibrate gamma
        radioButtons("decision_mode", "Threshold specification:",
                     choices = c("Specify posterior probability threshold γ" = "gamma",
                                 "Specify target Type-I error α" = "alpha"),
                     selected = "gamma"),

        conditionalPanel(
          condition = "input.decision_mode == 'gamma'",
          sliderInput(
            "gamma",
            label = HTML("Posterior probability threshold (\\( \\gamma \\))"),
            min = 80, max = 99, value = 90, step = 1, post = "%"
          ),
          uiOutput("decision_rule")   # placeholder for dynamic help text
        ),

        conditionalPanel(
          condition = "input.decision_mode == 'alpha'",
          sliderInput(
            "alpha",
            label = HTML("Target Type-I error (\\( \\alpha \\))"),
            value = 10, min = 1, max = 20, step = 1, post = "%"
          )
          ,
          helpText("γ will be calibrated so that the empirical Type-I error is ≈ α at the least-favourable null.")
        ),
      ),



      # --- Posterior Evaluation ---
      wellPanel(
        h4("Simulation settings"),
        numericInput("B", "Number of simulated trials (for Type-I & Power)", value = 2000, min = 100, step = 100),
        numericInput("ndraws", "Posterior draws per trial (for Pr(NI))", value = 2000, min = 200, step = 100),
        numericInput("seed", "Seed (optional)", value = 123, min = 1, step = 1),
        #helpText(HTML("Larger values give more precise results but increase runtime.")),
      ),

      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      h4("Operating characteristics at calibrated threshold"),
      textOutput("oc_text"),
      tableOutput("oc_table"),
      plotOutput("prni_plot", height = 350),
      helpText("Vertical line = calibrated γ. Shaded bars = distribution of posterior Pr(NI) across trials.",
               "Left panel: least-favourable null (Type-I). Right panel: assumed truth (Power).")
    )
  )
)

server <- function(input, output, session) {

  output$decision_rule <- renderUI({
    withMathJax(  # <- ensure newly injected HTML is typeset
      helpText(HTML(sprintf(
        "Trial is declared successful if the posterior probability that the treatment–control difference exceeds \\( \\Delta = %d\\%% \\) is at least \\( \\gamma = %d\\%% \\).",
        input$M, input$gamma
      )))
    )
  })

  # prior args as reactive list
  prior_args <- reactive({
    if (input$prior == "flat") return(list())
    list(
      a0 = input$a0 / 100,
      y_0 = input$y0,
      n_0 = input$n0,
      a_base = input$abase,
      b_base = input$bbase
    )
  })



  sim <- eventReactive(input$run, {
    pc <- input$pc / 100
    pt <- input$pt / 100
    M  <- input$M  / 100

    # Determine gamma depending on decision_mode
    if (input$decision_mode == "alpha") {
      alpha_target <- input$alpha / 100

      cal <- withProgress(
        message = "Calibrating γ...",
        detail  = "Running bisection search",
        value   = 0,
        {
          # You can simulate progress: each iteration increments progress
          # In bcts_calibrate_betaBinom_conj, add an argument like `progress_fun`
          # that calls `incProgress(1/maxit)` at each iteration.

          bcts_calibrate_betaBinom_conj(
            alpha   = alpha_target,
            p_c     = pc,  M = M,
            n_c     = input$nc,
            n_t     = input$nt,
            prior   = input$prior,
            prior_args = prior_args(),
            B_cal   = input$B,
            n_draws = input$ndraws,
            seed    = input$seed,
            show_progress = FALSE,   # suppress console bar
            verbose = FALSE,
            progress_fun = function(iter, maxit) {
              # iter=0 sent once before loop; avoid dividing by zero
              if (maxit > 0) incProgress(1 / maxit, detail = sprintf("Iteration %d of %d", max(1, iter), maxit))
            }
          )
        }
      )

      gamma_used <- cal$gamma
      cal_obj    <- cal
    } else {
      alpha_target <- NA_real_
      gamma_used <- input$gamma / 100
      cal_obj    <- NULL
    }

    # Type-I at LFN
    t1 <- bcts_type1_betaBinom_conj(
      B = input$B, p_c = pc, M = M,
      n_c = input$nc, n_t = input$nt,
      threshold = gamma_used,
      prior = input$prior, prior_args = prior_args(),
      n_draws = input$ndraws, show_progress = TRUE
    )

    # Power at assumed truth
    pw <- bcts_power_betaBinom_conj(
      B = input$B, p_c = pc, p_t = pt,
      n_c = input$nc, n_t = input$nt, M = M,
      threshold = gamma_used,
      prior = input$prior, prior_args = prior_args(),
      n_draws = input$ndraws, seed = input$seed, show_progress = FALSE
    )

    list(
      decision_mode = input$decision_mode,
      alpha_target  = alpha_target,
      gamma_used    = gamma_used,
      cal           = cal_obj,
      t1            = t1,   # keep full list (has CI)
      pw            = pw    # keep full list (has CI)
    )
  }, ignoreInit = TRUE)

  # small formatters
  fmt_pct <- function(x, d = 1) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("—")
    sprintf(paste0("%.", d, "f%%"), 100 * x)
  }
  fmt_ci <- function(lo, hi, d = 1) {
    if (is.null(lo) || is.null(hi) || length(lo) == 0 || length(hi) == 0 ||
        is.na(lo) || is.na(hi)) return("—")
    sprintf(paste0("[%.", d, "f, %.", d, "f]%%"), 100 * lo, 100 * hi)
  }

  output$oc_text <- renderText({
    s <- sim()
    if (is.null(s)) return("Run a simulation to see results.")

    # Build lines of text
    # Build lines of text
    lines <- c(
      sprintf("Gamma used: %.3f", s$gamma_used),
      sprintf(
        "Type-I error: %.1f%%  [%.1f%%, %.1f%%]  (MC SE ≈ %.1f%%)",
        100 * s$t1$estimate,
        100 * s$t1$ci_lower,
        100 * s$t1$ci_upper,
        100 * s$t1$mc_se
      ),
      sprintf("Power: %.1f%%  [%.1f%%, %.1f%%]",
              100 * s$pw$estimate,
              100 * s$pw$ci_lower,
              100 * s$pw$ci_upper)
    )

    paste(lines, collapse = "\n")
  })


}

shinyApp(ui, server)
