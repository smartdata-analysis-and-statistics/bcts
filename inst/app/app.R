# inst/app/app.R
# Shiny entrypoint for the bcts package app

# ---- Load packages that the APP needs ----
library(shiny)
library(ggplot2)
library(bcts)   # <- your package is built & installed during deploy

# ---- Locate the installed app directory (inside the bcts package) ----
app_dir <- system.file("app", package = "bcts")
stopifnot(nzchar(app_dir), dir.exists(app_dir))

# Optional: serve static assets from inst/app/www as /www
www_dir <- file.path(app_dir, "www")
if (dir.exists(www_dir)) addResourcePath("www", www_dir)

# Helper to source R files from inst/app/R and inst/app/modules
source_dir <- function(path) {
  if (dir.exists(path)) {
    fs <- list.files(path, pattern = "\\.R$", full.names = TRUE)
    for (f in fs) source(f, local = TRUE)
  }
}
source_dir(file.path(app_dir, "R"))
source_dir(file.path(app_dir, "modules"))


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
                     value = 2500, min = 100, step = 100),
        helpText(textOutput("design_mc_precision_text")),

        numericInput("ndraws", "Posterior draws per trial (for Pr(NI))", value = 2000, min = 200, step = 100),
        helpText(textOutput("post_mc_precision_text")),

        numericInput("seed", "Seed (optional)", value = 123, min = 1, step = 1),
        #helpText(HTML("Larger values give more precise results but increase runtime.")),
      ),



      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      mod_designsummary_ui("dsum"),
      mod_armpriors_ui("armpriors", height = 320),
      mod_oc_ui("oc"),

      # --- Sensitivity analysis ---
      # Sidebar
      mod_sensitivity_sidebar_ui("sens"),

      # Main panel (where you want outputs)
      mod_sensitivity_main_ui("sens", plot_height = 300)

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
  output$design_mc_precision_text <- renderText({
    B  <- as.numeric(req(input$B))
    se_max <- 0.5 / sqrt(B)                            # worst-case
    if (identical(input$decision_mode, "alpha")) {
      p <- input$alpha / 100
      se_alpha <- sqrt(p * (1 - p) / B)                # at alpha
      sprintf("Design MC SE: worst-case %.2f%% • at α=%d%%: %.2f%%",
              100*se_max, input$alpha, 100*se_alpha)
    } else {
      sprintf("Design MC SE (worst-case): %.2f%%", 100*se_max)
    }
  })

  output$post_mc_precision_text <- renderText({
    D <- as.numeric(req(input$ndraws))
    se_max <- 0.5 / sqrt(D)            # worst-case at p = 0.5
    hw_max <- 1.96 * se_max            # ~95% half-width
    sprintf("Posterior MC SE for Pr(NI) per trial: worst-case %.2f%% (±%.2f%%)",
            100*se_max, 100*hw_max)
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

  # results module consumes sim_res()
  mod_oc_server(
    id = "oc",
    sim = sim,                       # your eventReactive returning the results list
    decision_mode = reactive(input$decision_mode)
  )


  mod_designsummary_server(
    "dsum",
    pt = reactive(input$pt / 100),  nt = reactive(input$nt),
    pc = reactive(input$pc / 100),  nc = reactive(input$nc),
    M  = reactive(input$M  / 100),
    prior = reactive(input$prior),
    prior_args = reactive(list(
      a0 = if (input$prior == "power") input$a0/100 else 0,
      y_0 = if (input$prior == "power") input$y0 else 0,
      n_0 = if (input$prior == "power") input$n0 else 0,
      a_base = input$abase %||% 1,
      b_base = input$bbase %||% 1
    )),
    decision_mode = reactive(input$decision_mode),
    gamma = reactive(input$gamma / 100),
    alpha = reactive(input$alpha / 100),
    calibrate_on = reactive(input$calibrate_on),
    B = reactive(input$B),
    ndraws = reactive(input$ndraws),
    seed = reactive(input$seed)
  )



  mod_armpriors_server(
    "armpriors",
    p_t = reactive(input$pt / 100),
    n_t = reactive(input$nt),
    p_c = reactive(input$pc / 100),
    n_c = reactive(input$nc),
    prior = reactive(input$prior),
    prior_args = reactive(list(
      a_base = input$abase,
      b_base = input$bbase,
      a0     = if (input$prior == "power") input$a0 / 100 else 0,
      y_0    = if (input$prior == "power") input$y0 else 0,
      n_0    = if (input$prior == "power") input$n0 else 0
    ))
  )

  mod_sensitivity_server(
    "sens",
    sim        = sim,                           # your eventReactive from primary analysis
    pt         = reactive(input$pt / 100),
    nt         = reactive(input$nt),
    nc         = reactive(input$nc),
    M          = reactive(input$M / 100),
    prior      = reactive(input$prior),
    prior_args = reactive(prior_args()),
    B          = reactive(input$B),
    ndraws     = reactive(input$ndraws),
    seed       = reactive(input$seed),
    pc_current = reactive(input$pc)            # for pre-filling the range nicely
  )





}

shinyApp(ui, server)
