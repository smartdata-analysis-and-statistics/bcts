# modules/mod_designsummary.R
message("Loaded mod_designsummary.R")

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
        paste0(
          "<p><strong>Trial-level decision rule:</strong></p>",
          "<ul>",
          "<li>For a given simulated dataset, ",
          "\\(Pr(\\theta_t - \\theta_c > ", mjax_pct(vals$M, 0), ")\\) ",
          "is calculated from ", fmt_int(vals$ndraws), " posterior draws.</li>",
          "<li>The null hypothesis \\(H_0\\) is rejected if this probability ",
          "exceeds \\(\\gamma = ", mjax_pct(vals$gamma, 0), "\\).</li>",
          "</ul>"
        )
      } else {
        target <- switch(
          vals$calibrate_on,
          point = "the MC point estimate of Type-I error",
          upper = "the upper bound of the 95\\% MC CI for the Type-I error",
          lower = "the lower bound of the 95\\% MC CI for the Type-I error"
        )
        paste0(
          "<p><strong>Trial-level decision rule:</strong></p>",
          "<ul>",
          "<li>The posterior threshold \\(\\gamma\\) is calibrated so that ",
          target, " \\(\\approx \\alpha = ", mjax_pct(vals$alpha, 0), "\\).</li>",
          "<li>For a given simulated dataset, probabilities are calculated ",
          "from ", fmt_int(vals$ndraws), " posterior draws.</li>",
          "</ul>"
        )
      }

      sim_txt <- paste0(
        "<p><strong>Simulation setup:</strong></p>",
        "<ul>",
        "<li><b>Type-I error:</b> The trial-level decision rule is repeated ",
        fmt_int(vals$B), " times using datasets simulated under the least-favourable null.</li>",
        "<li><b>Power:</b> The same procedure is repeated ",
        fmt_int(vals$B), " times using datasets simulated under the assumed truths.</li>",
        "<li>Random seed = ", vals$seed, ".</li>",
        "</ul>"
      )

      # ---------------- NEW: Likelihood & Prior with configured numbers ----------------
      n_t <- vals$nt; n_c <- vals$nc

      dgm_eq <- paste0(
        "\\[\\begin{aligned}\n",
        "y_t &\\sim \\mathrm{Binomial}(", n_t, ",\\,", fmt_dec(vals$pt, 2), ") \\\\\n",
        "y_c &\\sim \\mathrm{Binomial}(", n_c, ",\\,", fmt_dec(vals$pc, 2), ")\n",
        "\\end{aligned}\\]"
      )


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
        dgm_eq,

        "<h5>Analysis model (likelihood & prior)</h5>",
        lik_eq, prior_eq,

        decision_txt,

        sim_txt,
        "</div>",
        "</details>"
      )

      withMathJax(HTML(html))
    })
  })
}
