# modules/mod_sensitivity.R

# ---- UI: sidebar controls ----
mod_sensitivity_sidebar_ui <- function(id) {
  ns <- NS(id)
  conditionalPanel(
    condition = sprintf("output['%s'] === true", ns("ready")),
    wellPanel(
      h4("Sensitivity analysis"),

      selectInput(ns("vary_which"), "Vary true response rate in:",
                  choices = c("Control group" = "control",
                              "Treatment group" = "treatment"),
                  selected = "control"
      ),

      conditionalPanel(
        condition = sprintf("input['%s'] === 'control'", ns("vary_which")),
        sliderInput(ns("theta_c_range"), "Control response rate (\\(\\theta_c\\))",
                    min = 0, max = 100, value = c(70, 95), step = 1, post = "%")
      ),

      conditionalPanel(
        condition = sprintf("input['%s'] === 'treatment'", ns("vary_which")),
        sliderInput(ns("theta_t_range"), "Treatment response rate (\\(\\theta_t\\))",
                    min = 0, max = 100, value = c(70, 95), step = 1, post = "%")
      ),

      numericInput(ns("grid"), "Number of grid points", value = 10, min = 3, step = 1, max = 15),
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
      condition = sprintf("output['%s'] === true", ns("has_results")),
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
                                   pt, pc, nt, nc,
                                   M,
                                   prior, prior_args,
                                   B, ndraws, seed,
                                   pc_current) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(x, y) if (is.null(x)) y else x
    .clamp <- function(x, a, b) pmin(pmax(x, a), b)

    observeEvent(sim(), {
      cat("Simulation changed — resetting sensitivity results\n")
      sens_stream(NULL)
    }, ignoreInit = TRUE)

    # Toggle visibility
    observe({
      toggle("theta_c_range", condition = input$vary_which == "control")
      toggle("theta_t_range", condition = input$vary_which == "treatment")
    })

    # Prefill both sliders to ensure input$theta_c_range and input$theta_t_range exist
    observe({
      updateSliderInput(session, "theta_c_range", value = c(75, 85))
      updateSliderInput(session, "theta_t_range", value = c(75, 85))
    })

    normalize_cols <- function(df) {
      nm <- names(df)
      nm <- gsub("\\.", "_", nm)   # just in case dots sneak in
      nm[nm == "lower"] <- "power_lo"
      nm[nm == "upper"] <- "power_hi"
      names(df) <- nm
      df
    }

    # ready flag
    ready <- reactive({ isTRUE(!is.null(sim())) })
    output$ready <- ready
    outputOptions(output, "ready", suspendWhenHidden = FALSE)

    # single definition (keep only this one)
    sens_stream <- reactiveVal(NULL)

    # results-present flag
    output$has_results <- reactive({
      df <- sens_stream()
      isTRUE(!is.null(df) && nrow(df) > 0)
    })
    outputOptions(output, "has_results", suspendWhenHidden = FALSE)

    # prefill θc range once when ready becomes TRUE (or is already TRUE at init)
    .prefilled <- reactiveVal(FALSE)

    # Prefill θc range (control)
    observeEvent(ready(), {
      req(ready(), !.prefilled())     # ensure sim() exists and run only once
      .prefilled(TRUE)

      s <- req(sim())

      # --- Control arm assumptions
      pc_assumed <- s$t1$settings$p_c %||% s$inputs$p_c %||% s$pc_assumed %||%
        (pc_current() / 100) %||% 0.85
      pc_assumed <- suppressWarnings(as.numeric(pc_assumed))
      if (!is.finite(pc_assumed)) pc_assumed <- 0.85
      pc_assumed <- .clamp(pc_assumed, 0, 1)

      pc_pct <- round(100 * pc_assumed)
      if (pc_pct == 0) { lo <- 0; hi <- 10 } else { lo <- max(0, pc_pct - 10); hi <- pc_pct }

      # update after panel renders
      session$onFlushed(function() {
        updateSliderInput(session, "theta_c_range", value = c(lo, hi))
      }, once = TRUE)
    }, ignoreInit = FALSE)

    # Prefill θt range (treatment)
    observeEvent(ready(), {
      req(ready())

      s <- req(sim())

      pt_assumed <- s$t1$settings$p_t %||% s$inputs$p_t %||% s$pt_assumed %||%
        (pt_current() / 100) %||% 0.85
      pt_assumed <- suppressWarnings(as.numeric(pt_assumed))
      if (!is.finite(pt_assumed)) pt_assumed <- 0.85
      pt_assumed <- .clamp(pt_assumed, 0, 1)

      pt_pct <- round(100 * pt_assumed)
      if (pt_pct == 0) {
        t_lo <- 0; t_hi <- 10
      } else {
        t_lo <- max(0, pt_pct - 10)
        t_hi <- pt_pct
      }

      session$onFlushed(function() {
        updateSliderInput(session, "theta_t_range", value = c(t_lo, t_hi))
      }, once = TRUE)
    }, ignoreInit = FALSE)

    # Run sensitivity with progress and streaming updates
    observeEvent(input$run, {
      s  <- req(sim())
      gam <- s$gamma_used
      vary_which <- input$vary_which %||% "control"

      cat("Triggered sensitivity run | vary_which =", vary_which, "\n")

      grid_n <- input$grid %||% 15
      Bv  <- B(); Nd <- ndraws(); Sd <- seed()
      ntv <- nt(); ncv <- nc()
      Mv <- M(); prv <- prior(); pav <- prior_args()

      # Determine min/max depending on what is varied
      if (vary_which == "control") {
        range_input <- req(input$theta_c_range)
        fixed_val <- pt()
        is_control <- TRUE
      } else {
        range_input <- req(input$theta_t_range)
        fixed_val <- pc()
        is_control <- FALSE
      }

      # Clamp to [0,1]
      min_val <- .clamp(min(range_input) / 100, 0, 1)
      max_val <- .clamp(max(range_input) / 100, 0, 1)
      if (max_val <= min_val) max_val <- min_val + 0.01

      pcs <- seq(min_val, max_val, length.out = grid_n)

      # What metrics to compute
      do_t1 <- "type1" %in% input$metrics
      do_pw <- "power" %in% input$metrics

      # Reset result stream
      sens_stream(NULL)

      # -- Before loop --
      cat("do_t1 =", do_t1, "| do_pw =", do_pw, "\n")
      cat("range_input =", paste(range_input, collapse = ", "), "\n")
      cat("pcs =", paste(round(pcs, 3), collapse = ", "), "\n")
      cat("Starting loop over", length(pcs), "points\n")

      withProgress(message = "Running sensitivity...", value = 0, {
        n <- length(pcs)
        for (i in seq_along(pcs)) {
          if (is_control) {
            pc <- pcs[i]
            pt <- fixed_val
          } else {
            pc <- fixed_val
            pt <- pcs[i]
          }

          t1 <- if (do_t1) bcts_type1_betaBinom_conj(
            B = Bv, p_c = pc, M = Mv, n_c = ncv, n_t = ntv,
            threshold = gam, prior = prv, prior_args = pav,
            n_draws = Nd, show_progress = FALSE
          ) else NULL

          pw <- if (do_pw) rct_beta_power(
            B = Bv, p_c = pc, p_t = pt, M = Mv, n_c = ncv, n_t = ntv,
            threshold = gam, prior = prv, prior_args = pav,
            n_draws = Nd, seed = Sd, show_progress = FALSE,
            method = "cpp"
          ) else NULL

          row <- data.frame(
            pc        = pc,
            pt        = pt,
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

          # Progress message
          rate_label <- if (is_control) "θc" else "θt"
          incProgress(1 / length(pcs),
                      detail = sprintf("Point %d/%d (%s = %.0f%%)",
                                       i, length(pcs), rate_label, 100 * pcs[i]))
        }
      })
    }, ignoreInit = TRUE)

    # --- Live plot (updates as sens_stream grows) ---------------------------
    output$plot <- renderPlot({
      df <- sens_stream()
      req(!is.null(df), nrow(df) > 0)
      s  <- req(sim())
      vary_which <- input$vary_which %||% "control"  # identify x-axis

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

      xvar <- if (vary_which == "control") df$pc else df$pt
      xlabel <- if (vary_which == "control") expression(theta[c]~"(true control rate)")
      else expression(theta[t]~"(true treatment rate)")


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
      long$x <- rep(xvar, each = 1)  # repeat per row (or match the length of each part)

      ggplot2::ggplot(long, ggplot2::aes(x = x, y = est)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 1.6, alpha = 0.85) +
        ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
        ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::labs(
          x = xlabel,
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


      out <- data.frame(
        `θ_c (true)` = sprintf("%.0f%%", 100*df$pc),
        `θ_t (true)` = sprintf("%.0f%%", 100*df$pt),
        check.names = FALSE, stringsAsFactors = FALSE
      )

      # Type-I error (combined column)
      if ("type1" %in% names(df) && any(!is.na(df$type1))) {
        out$`Type-I` <- ifelse(
          is.na(df$type1),
          "",
          sprintf("%.1f%% [%.1f%%, %.1f%%]",
                  100 * df$type1,
                  100 * df$type1_lo,
                  100 * df$type1_hi)
        )
      }

      # Power (combined column)
      if ("power" %in% names(df) && any(!is.na(df$power))) {
        out$`Power` <- ifelse(
          is.na(df$power),
          "",
          sprintf("%.1f%% [%.1f%%, %.1f%%]",
                  100 * df$power,
                  100 * df$power_lo,
                  100 * df$power_hi)
        )
      }
      out
    })
  })
}

