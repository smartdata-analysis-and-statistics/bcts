# modules/mod_designsummary.R

mod_designsummary_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      width = 12,
      wellPanel(
        htmlOutput(ns("design_text"))
      )
    )
  )
}

mod_designsummary_server <- function(
    id,
    # DGM (truths & sizes)
    pt, nt,
    pc, nc,
    # Analysis model
    prior,
    prior_args,              # list(a0, y_0, n_0, a_base, b_base)
    # Decision/calibration
    M, decision_mode, gamma, alpha, calibrate_on,
    # Simulation controls
    B, ndraws, seed
) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(x, y) if (is.null(x)) y else x
    fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
    fmt_int <- function(x) formatC(as.integer(x), big.mark = ",", format = "d")
    fmt_dec <- function(x, d = 3) formatC(x, format = "f", digits = d)
    mjax_pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "\\%")

    output$design_text <- renderUI({
      vals <- list(
        pt = pt(), nt = nt(),
        pc = pc(), nc = nc(),
        M  = M(),
        prior = prior(),
        pa = prior_args(),
        mode = decision_mode(),
        gamma = gamma(),
        alpha = alpha(),
        calibrate_on = calibrate_on(),
        B = B(), ndraws = ndraws(), seed = seed()
      )

      # Narrative bits (unchanged from your version)
      prior_txt <- if (identical(vals$prior, "power")) {
        sprintf(
          paste0(
            "Power prior on control: ",
            "<em>y</em><sub>0</sub> = %d of <em>n</em><sub>0</sub> = %d, ",
            "discount a<sub>0</sub> = %d%%, baseline Beta(%.1f, %.1f)."
          ),
          vals$pa$y_0 %||% 0, vals$pa$n_0 %||% 0, round(100 * (vals$pa$a0 %||% 0)),
          vals$pa$a_base %||% 1, vals$pa$b_base %||% 1
        )
      } else {
        sprintf("Baseline prior on control: Beta(%.1f, %.1f). Treatment prior: Beta(1,1).",
                vals$pa$a_base %||% 1, vals$pa$b_base %||% 1)
      }

      decision_txt <- if (identical(vals$mode, "gamma")) {
        sprintf(
          "Success if \\(\\Pr({\\theta}_t - {\\theta}_c > %s)\\ \\ge\\ %s\\).",
          mjax_pct(vals$M, 0), mjax_pct(vals$gamma, 0)
        )
      } else {
        target <- switch(vals$calibrate_on,
                         point = "the MC point estimate of Type-I error",
                         upper = "the upper bound of the 95\\% MC CI for the Type-I error",
                         lower = "the lower bound of the 95\\% MC CI for the Type-I error"
        )
        sprintf(
          "Posterior threshold \\(\\gamma\\) will be calibrated so that %s \\(\\approx\\) \\(\\alpha = %s\\) at the least-favourable null.",
          target, mjax_pct(vals$alpha, 0)
        )
      }

      sim_txt <- sprintf(
        "Planned Monte Carlo: %s trials (B) and %s posterior draws per trial; seed = %s.",
        fmt_int(vals$B), fmt_int(vals$ndraws), vals$seed
      )

      # ---------------- NEW: Likelihood & Prior with configured numbers ----------------
      n_t <- vals$nt; n_c <- vals$nc
      lik_eq <- paste0(
        "\\[\\begin{aligned}\n",
        "y_t\\mid\\theta_t &\\sim \\mathrm{Binomial}(", n_t, ",\\,\\theta_t)\\\\\n",
        "y_c\\mid\\theta_c &\\sim \\mathrm{Binomial}(", n_c, ",\\,\\theta_c)\n",
        "\\end{aligned}\\]"
      )

      abase <- vals$pa$a_base %||% 1
      bbase <- vals$pa$b_base %||% 1

      if (identical(vals$prior, "power")) {
        a0 <- vals$pa$a0   %||% 0
        y0 <- vals$pa$y_0  %||% 0
        n0 <- vals$pa$n_0  %||% 0
        a_pow <- abase + a0 * y0
        b_pow <- bbase + a0 * (n0 - y0)

        prior_eq <- paste0(
          "\\[\\begin{aligned}",

          # Treatment prior
          "\\theta_t &\\sim \\mathrm{Beta}(", fmt_dec(1, 3), ",\\,", fmt_dec(1, 3), ")\\\\\n",
          # Power prior as *numeric* Beta on control (after discounting history)
          "\\theta_c &\\sim \\mathrm{Beta}(",
          fmt_dec(a_pow, 3), ",\\,", fmt_dec(b_pow, 3), ")",
          "\\end{aligned}\\]"
        )
      } else {
        prior_eq <- paste0(
          "\\[\\begin{aligned}",
          "\\theta_t &\\sim \\mathrm{Beta}(", fmt_dec(1, 3), ",\\,", fmt_dec(1, 3), ")\\\\\n",
          "\\theta_c &\\sim \\mathrm{Beta}(", fmt_dec(abase, 3), ",\\,", fmt_dec(bbase, 3), ")",
          "\\end{aligned}\\]"
        )
      }
      # -------------------------------------------------------------------------------

      html <- paste0(
        "<details open><summary><h4>Design summary</h4></summary><div>",
        "<h5>Data-generating assumptions</h5>",
        "<ul>",
        "<li><b>Truths:</b> \\(\\theta_t = ", mjax_pct(vals$pt,0),
        "\\), \\(\\theta_c = ", mjax_pct(vals$pc,0), "\\)</li>",
        "<li><b>Sample sizes:</b> treatment \\(n_t = ", vals$nt,
        "\\), control \\(n_c = ", vals$nc, "\\)</li>",
        "</ul>",

        "<h5>Analysis model (likelihood & prior)</h5>",
        lik_eq, prior_eq,

        "<h5>Decision rule</h5><p>", decision_txt, "</p>",

        "<h5>Simulation setup</h5><p>", sim_txt, "</p>",
        "</div>",
        "</details>"
      )

      withMathJax(HTML(html))
    })
  })
}
