# modules/mod_sensitivity.R

# ---- UI: sidebar controls ----
mod_sensitivity_sidebar_ui <- function(id) {
  ns <- NS(id)
  wellPanel(
    h4("Sensitivity analysis"),

    # Message before primary analysis
    conditionalPanel(
      condition = sprintf("!output['%s']", ns("ready")),
      helpText("Run the primary analysis first to obtain γ; then you can run a sensitivity analysis.")
    ),

    # Controls shown only after primary analysis is ready
    conditionalPanel(
      condition = sprintf("output['%s']", ns("ready")),
      sliderInput(ns("pc_range"),
                  HTML("Vary control true rate (\\(\\theta_c\\))"),
                  min = 0, max = 100, value = c(70, 95), step = 1, post = "%"),
      numericInput(ns("grid"), "Number of grid points", value = 11, min = 3, step = 1),
      checkboxGroupInput(ns("metrics"), "Compute:",
                         choices  = c("Type-I at LFN" = "type1",
                                      "Power at assumed θt" = "power"),
                         selected = c("type1","power")),
      actionButton(ns("run"), "Run sensitivity", class = "btn btn-secondary")
    )
  )
}

# ---- UI: main-panel outputs ----
mod_sensitivity_main_ui <- function(id, plot_height = 300) {
  ns <- NS(id)
  tagList(
    conditionalPanel(
      condition = sprintf("output['%s']", ns("ready")),
      h4("Sensitivity results"),
      plotOutput(ns("plot"), height = plot_height),
      tableOutput(ns("table"))
    )
  )
}

# ---- Server ----
# sim: reactive returning the primary analysis object (must have $gamma_used)
# pt: reactive in [0,1] ; nt,nc: reactives (integers)
# M: reactive margin in [-1,1]
# prior: reactive "flat" or "power"
# prior_args: reactive list for your power prior settings
# B, ndraws, seed: reactives for MC settings
# pc_current: reactive current control slider (0-100) for prefill nicety
mod_sensitivity_server <- function(id,
                                   sim,
                                   pt, nt, nc,
                                   M,
                                   prior, prior_args,
                                   B, ndraws, seed,
                                   pc_current) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(x, y) if (is.null(x)) y else x

    normalize_cols <- function(df) {
      nm <- names(df)
      nm <- gsub("\\.", "_", nm)   # just in case dots sneak in
      nm[nm == "lower"] <- "power_lo"
      nm[nm == "upper"] <- "power_hi"
      names(df) <- nm
      df
    }

    # Gate UI until primary analysis exists
    output$ready <- reactive({ !is.null(sim()) })
    outputOptions(output, "ready", suspendWhenHidden = FALSE)

    # Prefill θc range around current slider when ready
    observeEvent(sim(), {
      pc <- round((pc_current() %||% 85))
      lo <- max(0, pc - 15); hi <- min(100, pc + 15)
      updateSliderInput(session, "pc_range", value = c(lo, hi))
    }, ignoreInit = TRUE)

    # Streaming results here
    sens_stream <- reactiveVal(NULL)

    # Run sensitivity with progress and streaming updates
    observeEvent(input$run, {
      s  <- req(sim())
      gam <- s$gamma_used

      # Build grid
      pc_min <- min(input$pc_range) / 100
      pc_max <- max(input$pc_range) / 100
      pcs    <- seq(pc_min, pc_max, length.out = input$grid)

      # Fixed knobs
      ntv <- nt(); ncv <- nc()
      ptv <- pt(); Mv <- M()
      prv <- prior(); pav <- prior_args()
      Bv  <- B(); Nd <- ndraws(); Sd <- seed()

      do_t1 <- "type1" %in% input$metrics
      do_pw <- "power" %in% input$metrics

      # reset stream
      sens_stream(NULL)

      withProgress(message = "Running sensitivity...", value = 0, {
        n <- length(pcs)
        for (i in seq_along(pcs)) {
          pc <- pcs[i]

          t1 <- if (do_t1) bcts_type1_betaBinom_conj(
            B = Bv, p_c = pc, M = Mv, n_c = ncv, n_t = ntv,
            threshold = gam, prior = prv, prior_args = pav,
            n_draws = Nd, show_progress = FALSE
          ) else NULL

          pw <- if (do_pw) bcts_power_betaBinom_conj(
            B = Bv, p_c = pc, p_t = ptv, M = Mv, n_c = ncv, n_t = ntv,
            threshold = gam, prior = prv, prior_args = pav,
            n_draws = Nd, seed = Sd, show_progress = FALSE
          ) else NULL

          row <- data.frame(
            pc        = pc,
            type1     = if (!is.null(t1)) t1$estimate else NA_real_,
            type1_lo  = if (!is.null(t1)) t1$ci_lower  else NA_real_,
            type1_hi  = if (!is.null(t1)) t1$ci_upper  else NA_real_,
            power     = if (!is.null(pw)) pw$estimate  else NA_real_,
            power_lo  = if (!is.null(pw)) pw$ci_lower  else NA_real_,
            power_hi  = if (!is.null(pw)) pw$ci_upper  else NA_real_,
            check.names = FALSE
          )
          row <- normalize_cols(row)     # <- force the names you expect

          cur <- sens_stream()
          cur <- if (is.null(cur)) row else rbind(cur, row)
          cur <- normalize_cols(cur)     # <- also normalize after binding
          sens_stream(cur)

          incProgress(1/n, detail = sprintf("Point %d/%d (θc = %.0f%%)", i, n, 100*pc))
        }
      })
    }, ignoreInit = TRUE)

    # --- Live plot (updates as sens_stream grows) ---------------------------
    output$plot <- renderPlot({
      df <- sens_stream()
      req(!is.null(df), nrow(df) > 0)
      s  <- req(sim())

      # helpers (keep if not already defined)
      get_col <- function(d, col, n) {
        if (col %in% names(d)) return(d[[col]])
        alt <- switch(col, power_lo = "lower", power_hi = "upper", col)
        if (alt %in% names(d)) return(d[[alt]])
        rep(NA_real_, n)
      }
      has_non_na <- function(x) any(!is.na(x))

      n <- nrow(df)

      # safely pull columns (NA vectors if missing)
      t1_est <- get_col(df, "type1",    n)
      t1_lo  <- get_col(df, "type1_lo", n)
      t1_hi  <- get_col(df, "type1_hi", n)
      pw_est <- get_col(df, "power",    n)
      pw_lo  <- get_col(df, "power_lo", n)
      pw_hi  <- get_col(df, "power_hi", n)

      parts <- list()
      if (has_non_na(t1_est)) {
        parts$type1 <- data.frame(
          pc = df$pc, metric = factor("Type-I", levels = c("Type-I","Power")),
          est = t1_est, lo = t1_lo, hi = t1_hi
        )
      }
      if (has_non_na(pw_est)) {
        parts$power <- data.frame(
          pc = df$pc, metric = factor("Power", levels = c("Type-I","Power")),
          est = pw_est, lo = pw_lo, hi = pw_hi
        )
      }
      validate(need(length(parts) > 0, "Select at least one metric to compute."))

      long <- do.call(rbind, parts)

      ggplot2::ggplot(long, ggplot2::aes(x = pc, y = est)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 1.6, alpha = 0.85) +
        ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
        ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::labs(
          x = expression(theta[c]~"(true control rate)"),
          y = "Operating characteristic",
          subtitle = sprintf("γ used = %.3f, Δ = %d%%, n_t = %d, n_c = %d, B = %d",
                             s$gamma_used, round(100*M()), nt(), nc(), B())
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          strip.text = ggplot2::element_text(face = "bold"),
          legend.position = "none"
        )
    })

    # --- Table (robust to deselected metrics) -------------------------------
    output$table <- renderTable({
      df <- sens_stream()
      req(!is.null(df), nrow(df) > 0)

      n <- nrow(df)
      out <- data.frame(`θ_c (true)` = sprintf("%.0f%%", 100*df$pc),
                        check.names = FALSE, stringsAsFactors = FALSE)

      if ("type1" %in% names(df) && any(!is.na(df$type1))) {
        out$`Type-I`        <- ifelse(is.na(df$type1), "", sprintf("%.1f%%", 100*df$type1))
        out$`Type-I 95% CI` <- ifelse(is.na(df$type1_lo), "",
                                      sprintf("[%.1f%%, %.1f%%]", 100*df$type1_lo, 100*df$type1_hi))
      }
      if ("power" %in% names(df) && any(!is.na(df$power))) {
        out$`Power`         <- ifelse(is.na(df$power), "", sprintf("%.1f%%", 100*df$power))
        out$`Power 95% CI`  <- ifelse(is.na(df$power_lo), "",
                                      sprintf("[%.1f%%, %.1f%%]", 100*df$power_lo, 100*df$power_hi))
      }
      out
    })
  })
}

# tiny helper (keeps module self-contained)
`%||%` <- function(x, y) if (is.null(x)) y else x
