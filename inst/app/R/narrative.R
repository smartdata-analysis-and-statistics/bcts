design_narrative <- function(input_vals, s) {
  if (is.null(s)) return(HTML("<p>Run a simulation to generate the design summary.</p>"))
  with(input_vals, {
    prior_txt <- if (prior == "flat") {
      "Flat prior on both arms (no external information)."
    } else {
      sprintf("Power prior: y0=%d, n0=%d, a0=%d%%, Beta(%.1f, %.1f).",
              y0, n0, a0, abase, bbase)
    }
    decision_txt <- if (mode == "gamma") {
      sprintf("Decision: Pr(θt − θc > Δ) ≥ γ with Δ=%s, γ=%s.",
              fmt_pct(M,0), fmt_pct(gamma,0))
    } else {
      sprintf("γ calibrated to α=%s (%s). Calibrated γ=%s.",
              fmt_pct(alpha,0),
              switch(calibrate_on, point="point", upper="upper 95% CI", lower="lower 95% CI"),
              fmt_prob(s$gamma_used,3))
    }
    HTML(paste0(
      "<ul>",
      "<li><b>Margin</b>: Δ=", fmt_pct(M,0), "</li>",
      "<li><b>Sample sizes</b>: nt=", nt, ", nc=", nc, "</li>",
      "<li><b>Prior</b>: ", prior_txt, "</li>",
      "<li><b>", decision_txt, "</b></li>",
      "</ul>"
    ))
  })
}
