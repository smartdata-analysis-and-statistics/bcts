# modules/mod_designsummary.R

mod_designsummary_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Design summary"),
    htmlOutput(ns("design_text"))
  )
}

mod_designsummary_server <- function(
    id,
    pt, nt,                      # reactives: treatment truth in [0,1], n_t
    pc, nc,                      # reactives: control truth in [0,1], n_c
    M,                           # reactive: margin (risk-diff) in [-1,1]
    prior,                       # reactive: "flat" or "power"
    prior_args,                  # reactive: list(a0 in [0,1], y_0, n_0, a_base, b_base)
    decision_mode,               # reactive: "gamma" or "alpha"
    gamma, alpha, calibrate_on,  # reactives: gamma/alpha in [0,1], calibrate_on: point/upper/lower
    B, ndraws, seed              # reactives: simulation knobs
) {
  moduleServer(id, function(input, output, session) {

    # tiny helper
    `%||%` <- function(x, y) if (is.null(x)) y else x
    fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
    fmt_int <- function(x) formatC(as.integer(x), big.mark = ",", format = "d")

    output$design_text <- renderUI({
      # pull current inputs (no dependency on sim())
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

      # prior paragraph
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

      # decision paragraph
      # helper: percent for MathJax
      mjax_pct <- function(x, d = 0) {
        # x in [0,1]; returns like "20\\%"
        paste0(formatC(100 * x, format = "f", digits = d), "\\%")
      }

      # inside your module/server:
      decision_txt <- if (identical(vals$mode, "gamma")) {
        sprintf(
          "Success if \\(\\Pr(\\hat{\\theta}_t - \\hat{\\theta}_c > %s)\\ \\ge\\ %s\\).",
          mjax_pct(vals$M, 0),     # e.g. "-20\\%"
          mjax_pct(vals$gamma, 0)  # e.g. "90\\%"
        )
      }  else {
        target <- switch(vals$calibrate_on,
                         point = "the MC point estimate of Type-I error",
                         upper = "the upper bound of the 95% MC CI for the Type-I error",
                         lower = "the lower bound of the 95% MC CI for the Type-I error")
        sprintf(
          "Posterior threshold \\(\\gamma\\) will be calibrated so that %s ≈ \\(\\alpha = %s\\) at the least-favourable null.",
          target, fmt_pct(vals$alpha, 0)
        )
      }

      # simulation knobs (informational only)
      sim_txt <- sprintf(
        "Planned Monte Carlo: %s trials (B) and %s posterior draws per trial; seed = %s.",
        fmt_int(vals$B), fmt_int(vals$ndraws), vals$seed
      )

      # assemble (design only; no results)
      html <- paste0(
        "<p><strong>Design:</strong> Binary endpoint with Beta–Binomial conjugate updates.</p>",
        "<ul>",
          "<li><b>Assumed truths</b>: ",
          "\\(\\theta_t = ", mjax_pct(vals$pt, 0), "\\), ",
          "\\(\\theta_c = ", mjax_pct(vals$pc, 0), "\\)",
          "</li>",
        "<li><b>Sample sizes</b>: treatment \\(n_t = ", vals$nt,
        "\\), control \\(n_c = ", vals$nc, "\\)</li>",
        "<li><b>Margin</b>: \\(\\Delta = ", mjax_pct(vals$M, 0), "\\) ",
        "(Δ &lt; 0 non-inferiority; Δ ≥ 0 superiority)</li>",
        "<li><b>Prior</b>: ", prior_txt, "</li>",
        "<li><b>Decision rule</b>: ", decision_txt, "</li>",
        "<li><b>Simulation setup</b>: ", sim_txt, "</li>",
        "</ul>"
      )

      withMathJax(HTML(html))
    })
  })
}
