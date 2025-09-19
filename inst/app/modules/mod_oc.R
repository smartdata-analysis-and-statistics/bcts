# modules/mod_oc.R
# Operating characteristics (table + optional calibration plot)

# modules/mod_oc.R
mod_oc_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("oc_box"))   # dynamically rendered box
}

mod_oc_server <- function(id, sim, decision_mode) {
  moduleServer(id, function(input, output, session) {

    # small formatters
    fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
    fmt_ci  <- function(lo, hi, d = 1) sprintf(paste0("[%.", d, "f%%, %.", d, "f%%]"),
                                               100 * lo, 100 * hi)

    # Whole OC panel (only shown if sim() has results)
    output$oc_box <- renderUI({
      s <- sim()
      if (is.null(s)) return(NULL)

      ns <- session$ns
      wellPanel(
        h4("Operating characteristics"),
        textOutput(ns("blurb")),
        tableOutput(ns("table")),
        conditionalPanel(
          condition = sprintf("output['%s']", ns("show_cal")),
          h4("Calibration: Type-I error across γ tried"),
          plotOutput(ns("cal_trace_plot"), height = 300)
        )
      )
    })

    # --- Blurb ----------------------------------------------------------------
    output$blurb <- renderText({
      mode <- req(decision_mode())
      if (identical(mode, "alpha")) {
        "Simulation results showing the calibrated posterior threshold (γ) and corresponding Type-I error at the least-favourable null, as well as the power at the assumed truth."
      } else {
        "Simulation results showing the specified posterior threshold (γ), with Type-I error at the least-favourable null and power at the assumed truth."
      }
    })

    # --- OC table --------------------------------------------------------------
    output$table <- renderTable({
      s <- sim()
      mode <- req(decision_mode())
      req(!is.null(s))

      # Determine if we should use calibration info
      use_cal <- identical(mode, "alpha") && !is.null(s$cal)

      # --- Type-I Error ---
      type1_est  <- if (use_cal) fmt_pct(s$cal$type1["estimate"], 1) else fmt_pct(s$t1$estimate, 1)
      type1_ci   <- if (use_cal) fmt_ci(s$cal$type1["ci_lower"], s$cal$type1["ci_upper"], 1) else fmt_ci(s$t1$ci_lower, s$t1$ci_upper, 1)
      type1_mcse <- if (use_cal) fmt_pct(s$cal$type1["mc_se"], 1) else fmt_pct(s$t1$mc_se, 1)

      # --- Power (always from s$pw) ---
      power_est   <- fmt_pct(s$pw$estimate, 1)
      power_ci    <- fmt_ci(s$pw$ci_lower, s$pw$ci_upper, 1)
      power_mcse  <- fmt_pct(s$pw$mc_se, 1)

      # --- Build table ---
      data.frame(
        Metric   = c("Gamma used", "Type-I error", "Power"),
        Estimate = c(sprintf("%.3f", s$gamma_used), type1_est, power_est),
        `95% CI` = c("—", type1_ci, power_ci),
        `MC SE`  = c("—", type1_mcse, power_mcse),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    # Expose whether calibration plot should be shown
    output$show_cal <- reactive({
      mode_is_alpha <- identical(decision_mode(), "alpha")
      s <- sim()
      mode_is_alpha && !is.null(s) && !is.null(s$cal)
    })
    outputOptions(output, "show_cal", suspendWhenHidden = FALSE)

    # Calibration trace plot
    output$cal_trace_plot <- renderPlot({
      s <- sim()
      validate(need(!is.null(s) && !is.null(s$cal),
                    "Calibration plot available only when γ was calibrated."))

      plot(s$cal)  # Calls the S3 method plot.bcts_calibration()
    })
  })
}
