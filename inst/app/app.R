# app.R
# Minimal Shiny app for Bayesian NI (Beta–Binomial, conjugate sampling)

library(shiny)
library(ggplot2)

# Source helpers & modules
lapply(list.files("R", full.names = TRUE), source)
lapply(list.files("modules", full.names = TRUE), source)


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
      tagList(
        # Panel 1: Control arm (always shown)
        wellPanel(
          h4("Control Arm Assumptions"),
          sliderInput("pc", HTML("True response rate (\\(\\theta_c\\))"),
                      min = 0, max = 100, value = 85, step = 1, post = "%"),
          numericInput("nc", HTML("Number of randomized patients (\\(n_c\\))"),
                       value = 29, min = 1, step = 1),

          selectInput(
            "prior",
            "Prior distribution",
            choices = c("Flat (no external evidence)" = "flat",
                        "Power prior (with historical data)" = "power"),
            selected = "flat"
          )
        ),

        # Panel 2: External data config (only when prior == power)
        conditionalPanel(
          condition = "input.prior == 'power'",
          wellPanel(
            h4("External Data (Power Prior)"),
            sliderInput("a0", "Discount factor a₀ (0 = ignore history, 1 = full borrow)",
                        min = 0, max = 100, value = 50, step = 1, post = "%"),
            numericInput("y0", "Historical responders (y₀)", value = 64, min = 0, step = 1),
            numericInput("n0", "Historical sample size (n₀)", value = 75, min = 1, step = 1),
            numericInput("abase", "Baseline Beta a_base", value = 1, min = 0.001, step = 0.1),
            numericInput("bbase", "Baseline Beta b_base", value = 1, min = 0.001, step = 0.1)
          )
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
          ),
          selectInput(
            "calibrate_on",
            "Calibrate Type-I on:",
            choices = c(
              "Point estimate Pr(reject | H₀)"      = "point",
              "Upper 95% MC CI (conservative)"      = "upper",
              "Lower 95% MC CI (liberal)"           = "lower"
            ),
            selected = "upper"
          ),
          helpText("γ will be calibrated so that the chosen Type-I metric ≈ α (within tolerance) at the least-favourable null.")
        )
      ),



      # --- Posterior Evaluation ---
      wellPanel(
        h4("Simulation settings"),
        numericInput("B",
                     "Number of simulated trials (for Type-I & Power)",
                     value = 2000, min = 100, step = 100),
        helpText(textOutput("mc_precision_text")),

        numericInput("ndraws", "Posterior draws per trial (for Pr(NI))", value = 2000, min = 200, step = 100),
        numericInput("seed", "Seed (optional)", value = 123, min = 1, step = 1),
        #helpText(HTML("Larger values give more precise results but increase runtime.")),
      ),

      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      h4("Design summary"),
      htmlOutput("design_text"),
      hr(),
      h4("Operating characteristics"),
      textOutput("oc_text"),
      tableOutput("oc_table"),

      # --- Calibration plot only when decision_mode == 'alpha' ---
      conditionalPanel(
        condition = "input.decision_mode == 'alpha'",
        h4("Calibration: Type-I error across γ tried"),
        plotOutput("cal_trace_plot", height = 300)
      )
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

  # Server (change renderUI -> renderText)
  output$mc_precision_text <- renderText({
    B  <- as.numeric(req(input$B))
    se <- 0.5 / sqrt(B)   # worst-case MC SE
    sprintf("Max MC SE: %.2f%%", 100 * se)
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
    # scale to probabilities
    pc <- input$pc / 100
    pt <- input$pt / 100
    M  <- input$M  / 100

    withProgress(message = "Calibrating / simulating…", value = 0, {
      # Map bisection progress to Shiny progress bar (only used if mode == "alpha")
      pf <- function(iter, maxit) {
        if (maxit <= 0) return()
        frac <- max(0, min(1, iter / maxit))
        setProgress(frac, detail = sprintf("Calibration %d / %d", max(1L, iter), maxit))
      }

      out <- run_all_oc(
        pt = pt, nt = input$nt,
        pc = pc, nc = input$nc,
        M  = M,
        mode         = input$decision_mode,          # "gamma" or "alpha"
        gamma        = input$gamma / 100,            # ignored if mode == "alpha"
        alpha        = input$alpha / 100,            # used only if mode == "alpha"
        calibrate_on = input$calibrate_on,
        prior        = input$prior,
        prior_args   = prior_args(),
        B            = input$B,                      # your UI uses id "B"
        ndraws       = input$ndraws,
        seed         = input$seed,
        progress_fun = pf,                           # will be used by run_all_oc()
        verbose      = FALSE
      )

      setProgress(1, detail = "Done")
      out
    })
  }, ignoreInit = TRUE)

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

  output$cal_trace_plot <- renderPlot({
    s <- sim()
    # only show when we actually calibrated
    if (is.null(s) || is.null(s$cal)) return(NULL)

    tr <- s$cal$trace            # data.frame with iter, gamma_try, type1, lo, hi, diff
    alpha <- s$cal$alpha
    gamma <- s$gamma_used

    # basic sanity check
    if (!all(c("gamma_try","type1") %in% names(tr))) return(NULL)

    ggplot(tr, aes(x = gamma_try, y = type1)) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.002, alpha = 0.4) +
      geom_point(size = 2) +
      geom_line() +
      geom_hline(yintercept = alpha, linetype = "dotted", color = "red") +
      geom_vline(xintercept = gamma, linetype = "dashed", color = "blue") +
      labs(
        x = expression(gamma),
        y = "Estimated Type-I error",
        subtitle = "Points are bisection iterations with 95% CI; dashed = calibrated γ, dotted = target α"
      ) +
      theme_minimal(base_size = 12)
  })

  output$oc_table <- renderTable({
    s <- sim()
    if (is.null(s)) return(NULL)

    df <- data.frame(
      Metric = c("Gamma used", "Type-I error", "Power"),
      Estimate = c(
        sprintf("%.3f", s$gamma_used),
        sprintf("%.1f%%", 100 * s$t1$estimate),
        sprintf("%.1f%%", 100 * s$pw$estimate)
      ),
      `95% CI` = c(
        "—",
        sprintf("[%.1f%%, %.1f%%]", 100 * s$t1$ci_lower, 100 * s$t1$ci_upper),
        sprintf("[%.1f%%, %.1f%%]", 100 * s$pw$ci_lower, 100 * s$pw$ci_upper)
      ),
      `MC SE` = c(
        "—",
        sprintf("%.1f%%", 100 * s$t1$mc_se),
        "—"
      ),
      check.names = FALSE,  # <- keep column labels as written
      stringsAsFactors = FALSE
    )

    df
  })

  output$design_text <- renderUI({
    s <- NULL
    if (!is.null(isolate(reactiveValuesToList(input)$run)) && !is.null(sim())) s <- sim()

    input_vals <- list(
      prior        = input$prior,
      y0           = input$y0,
      n0           = input$n0,
      a0           = input$a0,
      abase        = input$abase,
      bbase        = input$bbase,
      mode         = input$decision_mode,
      gamma        = input$gamma / 100,
      alpha        = input$alpha / 100,
      calibrate_on = input$calibrate_on,
      nt           = input$nt,
      nc           = input$nc,
      M            = input$M / 100
    )

    withMathJax(design_narrative(input_vals, s))
  })





}

shinyApp(ui, server)
