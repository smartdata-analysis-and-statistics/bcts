# inst/app/app.R
# Shiny entrypoint for the bcts package app

# ---- Load packages that the APP needs ----
library(shiny)
library(ggplot2)
library(bcts)

# Helper to source R files from inst/app/R and inst/app/modules
source_dir <- function(path) {
  if (dir.exists(path)) {
    fs <- list.files(path, pattern = "\\.R$", full.names = TRUE)
    for (f in fs) {
      message("Sourcing: ", basename(f))
      tryCatch({
        source(f, local = globalenv())  # Use globalenv for visibility
      }, error = function(e) {
        message("Error sourcing ", f, ": ", e$message)
      })
    }
  }
}

USE_LOCAL <- interactive() && dir.exists(file.path("inst", "app"))  # Only TRUE when running locally

app_dir <- if (USE_LOCAL) {
  file.path("inst", "app")
} else {
  "."
}

source_dir(file.path(app_dir, "R"))
source_dir(file.path(app_dir, "modules"))


# If your functions live in a package, uncomment:
# library(bcts)

ui <- navbarPage(
  title = "Bayesian Trial Simulation (Beta–Binomial, conjugate)",
  tabPanel(
    title = "Randomized trial",
    fluidPage(
      withMathJax(),   # enable LaTeX rendering for the whole page

      # --- Description row at top (full width) ---
      fluidRow(
        column(
          width = 12,
          div(
            style = "margin-bottom: 20px;",
            h4("Simulate binary responder outcomes from a randomized trial with a treatment and control arm."),
            p("Use the controls in the left panel to define the data-generating assumptions, prior distributions, and decision criteria.
            The design summary and operating characteristics will be shown on the right.")
          )
        )
      ),


      sidebarPanel(

        wellPanel(
          h4("Summary Output Type"),
          radioButtons(
            inputId = "summary_type",
            label = "Choose summary format:",
            choices = c("Narrative" = "narrative", "Technical" = "technical"),
            selected = "narrative"
          )
        ),

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

        decision_criteria_ui(),

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



        actionButton("run", "Run simulation", class = "btn-primary"),

        hr(),
        tags$small(
          paste("bcts version:", utils::packageVersion("bcts"))
        )
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.summary_type == 'technical'",
          mod_designsummary_ui("dsum")
        ),
        conditionalPanel(
          condition = "input.summary_type == 'narrative'",
          mod_narrative_ui("narrative")
        ),
        mod_armpriors_ui("armpriors", height = 320),
        mod_oc_ui("oc"),

        # --- Sensitivity analysis ---
        # Sidebar
        mod_sensitivity_sidebar_ui("sens"),

        # Main panel (where you want outputs)
        mod_sensitivity_main_ui("sens", plot_height = 300)


      )
    )
  ),

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

          sliderInput("pt_sa", HTML("True response rate (\\( \\theta_t \\))"),
                      min = 0, max = 100, value = 80, step = 1, post = "%"),

          numericInput("nt_sa", HTML("Sample size (\\( n_t \\))"),
                       value = 40, min = 1, step = 1),

          selectInput(
            "prior_sa",
            "Prior distribution",
            choices = c("Flat (Beta(1,1))" = "flat",
                        "Custom Beta prior" = "beta"),
            selected = "flat"
          ),

          conditionalPanel(
            condition = "input.prior_sa == 'beta'",
            numericInput("abase_sa", "Prior a_base", value = 1, min = 0.01, step = 0.1),
            numericInput("bbase_sa", "Prior b_base", value = 1, min = 0.01, step = 0.1)
          ),

          sat_decision_criteria_ui("crit_sa"),
        ),

        mainPanel(
          verbatimTextOutput("sa_summary"),
          plotOutput("sa_power_plot")
        )
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
  output$`crit_sa-decision_rule` <- renderUI({
    withMathJax(
      helpText(HTML(sprintf(
        "Trial is declared successful if posterior \\( \\Pr(\\theta > %.2f) \\geq %.1f\\%% \\).",
        input[["crit_sa-M_sa"]] / 100,
        input[["crit_sa-gamma_sa"]]
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

  mod_narrative_server(
    "narrative",
    pt = reactive(input$pt / 100),  nt = reactive(input$nt),
    pc = reactive(input$pc / 100),  nc = reactive(input$nc),
    M  = reactive(input$M  / 100),
    prior = reactive(input$prior),
    prior_args = prior_args,
    decision_mode = reactive(input$decision_mode),
    gamma = reactive(input$gamma / 100),
    alpha = reactive(input$alpha / 100),
    B = reactive(input$B),
    ndraws = reactive(input$ndraws),
    seed = reactive(input$seed)
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

  sa_results <- reactive({
    req(input$pt_sa, input$nt_sa,
        input[["crit_sa-M_sa"]],
        input[["crit_sa-gamma_sa"]])

    prior_type <- input$prior_sa
    a_base <- if (prior_type == "beta") input$abase_sa else 1
    b_base <- if (prior_type == "beta") input$bbase_sa else 1

    B <- input$B_sa
    n_draws <- input$ndraws_sa
    seed <- input$seed_sa

    pt <- input$pt_sa / 100
    M <- input[["crit_sa-M_sa"]] / 100
    gamma <- input[["crit_sa-gamma_sa"]] / 100

    power_res <- bcts::singlearm_beta_power(
      B = B,
      p_t = pt,
      n_t = input$nt_sa,
      M = M,
      threshold = gamma,
      prior = prior_type,
      a_base = a_base,
      b_base = b_base,
      method = "exact",
      show_progress = FALSE
    )

    type1_res <- bcts::singlearm_beta_type1(
      B = B,
      n_t = input$nt_sa,
      M = M,
      threshold = gamma,
      prior = prior_type,
      a_base = a_base,
      b_base = b_base,
      n_draws = n_draws,
      method = "exact",
      show_progress = FALSE
    )

    list(power = power_res, type1 = type1_res)
  })

  output$sa_summary <- renderPrint({
    res <- sa_results()
    power_res <- res$power
    type1_res <- res$type1

    cat("POWER ANALYSIS\n")
    cat(sprintf("Estimated power: %.2f%%\n", 100 * power_res$estimate))
    cat(sprintf("MC standard error: %.2f%%\n", 100 * power_res$mc_se %||% NA_real_))
    cat(sprintf("Successes: %d out of %d simulations\n\n",
                power_res$successes %||% NA_integer_,
                power_res$B %||% NA_integer_))

    cat("TYPE-I ERROR ANALYSIS\n")
    cat(sprintf("Estimated Type-I error: %.2f%%\n", 100 * type1_res$estimate))
    cat(sprintf("MC standard error: %.2f%%\n", 100 * type1_res$mc_se %||% NA_real_))
    cat(sprintf("False positives: %d out of %d simulations\n",
                type1_res$successes %||% NA_integer_,
                type1_res$B %||% NA_integer_))
  })

}

tryCatch(
  {
    shiny::shinyApp(
      ui = ui,
      server = server
    )
  },
  error = function(e) {
    message("🚨 App failed to launch: ", conditionMessage(e))
    stop(e)
  }
)
