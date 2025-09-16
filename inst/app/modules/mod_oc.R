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
      req(!is.null(s))

      data.frame(
        Metric   = c("Gamma used", "Type-I error", "Power"),
        Estimate = c(
          sprintf("%.3f", s$gamma_used),
          fmt_pct(s$t1$estimate, 1),
          fmt_pct(s$pw$estimate, 1)
        ),
        `95% CI` = c(
          "—",
          fmt_ci(s$t1$ci_lower, s$t1$ci_upper, 1),
          fmt_ci(s$pw$ci_lower, s$pw$ci_upper, 1)
        ),
        `MC SE`  = c(
          "—",
          fmt_pct(s$t1$mc_se, 1),
          fmt_pct(s$pw$mc_se, 1)
        ),
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

      tr <- s$cal$trace
      alpha <- s$cal$alpha
      gamma <- s$gamma_used

      ggplot2::ggplot(tr, ggplot2::aes(x = gamma_try, y = type1)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                               width = 0.002, alpha = 0.4) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = alpha, linetype = "dotted", color = "red") +
        ggplot2::geom_vline(xintercept = gamma, linetype = "dashed", color = "blue") +
        ggplot2::labs(
          x = expression(gamma),
          y = "Estimated Type-I error",
          subtitle = "Points are bisection iterations with 95% CI; dashed = calibrated γ, dotted = target α"
        ) +
        ggplot2::theme_minimal(base_size = 12)
    })
  })
}
