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
          numericInput(
            "alpha",
            label = HTML("Target Type-I error (\\( \\alpha \\))"),
            value = 0.10, min = 0.001, max = 0.5, step = 0.001
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

      hr(),
      h4("Decision threshold & operating characteristics"),

      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      h4("Operating characteristics at calibrated threshold"),
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

  # Build prior args list (reactive)
  prior_args <- reactive({
    if (input$prior == "flat") return(list())
    list(a0 = input$a0/100, y_0 = input$y0, n_0 = input$n0, a_base = input$abase, b_base = input$bbase)
  })

  sim <- eventReactive(input$run, {
    pc <- input$pc / 100
    pt <- input$pt / 100


    # 2) Estimate Type-I and Power at calibrated gamma (same B)
    t1 <- bayesNI_type1_betaBinom_conj(
      B = input$B, p_c = pc, M = input$M/100,
      n_c = input$nc, n_t = input$nt,
      threshold = input$gamma/100,
      prior = input$prior, prior_args = prior_args(),
      n_draws = input$ndraws, show_progress = FALSE
    )

    pw <- bayesNI_power_betaBinom_conj(
      B = input$B, p_c = pc, p_t = pt,
      n_c = input$nc, n_t = input$nt, M = input$M/100,
      threshold =  input$gamma/100,
      prior = input$prior, prior_args = prior_args(),
      n_draws = input$ndraws, seed = input$seed, show_progress = FALSE
    )

    # 3) Distributions of posterior Pr(NI) for plotting
    pr_type1 <- prNI_draws_conj(
      B = input$B, p_c = pc, p_t = pc + input$M/100, n_c = input$nc, n_t = input$nt, M = input$M/100,
      prior = input$prior, prior_args = prior_args(),
      n_draws = input$ndraws, seed = input$seed
    )
    pr_power <- prNI_draws_conj(
      B = input$B, p_c = pc, p_t = pt, n_c = input$nc, n_t = input$nt, M = input$M/100,
      prior = input$prior, prior_args = prior_args(),
      n_draws = input$ndraws, seed = input$seed
    )

    list(gamma =  input$gamma/100, type1 = t1$type1, power = pw,
         pr_type1 = pr_type1, pr_power = pr_power)
  }, ignoreInit = TRUE)

  output$oc_table <- renderTable({
    s <- sim(); if (is.null(s)) return(NULL)
    data.frame(
      Metric   = c("Gamma", "Type-I error", "Power"),
      Estimate = c(
        sprintf("%.3f", s$gamma),
        sprintf("%.3f", s$type1),
        sprintf("%.3f", s$power)
      ),
      check.names = FALSE
    )
  })

  output$prni_plot <- renderPlot({
    s <- sim(); if (is.null(s)) return(NULL)
    df <- rbind(
      data.frame(pr = s$pr_type1, panel = "LFN (Type-I)"),
      data.frame(pr = s$pr_power, panel = "Truth (Power)")
    )
    ggplot(df, aes(pr)) +
      geom_histogram(bins = 40) +
      geom_vline(xintercept = s$cal$gamma, linetype = "dashed") +
      facet_wrap(~panel, nrow = 1) +
      labs(x = "Posterior Pr(NI)", y = "Count")
  })
}

shinyApp(ui, server)
