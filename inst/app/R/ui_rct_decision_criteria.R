# inst/app/R/ui_rct_decision_criteria.R

decision_criteria_ui <- function(id = NULL) {
  ns <- NS(id)  # allows modular usage, or just use identity

  wellPanel(
    h4("Decision criteria"),
    sliderInput(
      ns("M"),
      label = HTML("Decision Margin (\\( \\Delta \\))"),
      min   = -100, max = 100, value = -20, step = 1, post = "%"
    ),
    helpText(
      "Δ < 0: non-inferiority (treatment may be up to |Δ| worse).",
      "Δ ≥ 0: superiority (treatment must be at least Δ better).",
      "Assumes higher response rates are better (responder events)."
    ),

    radioButtons(ns("decision_mode"), "Threshold specification:",
                 choices = c(
                   "Specify posterior probability threshold γ" = "gamma",
                   "Specify target Type-I error α" = "alpha"
                 ),
                 selected = "gamma"),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'gamma'", ns("decision_mode")),
      sliderInput(
        ns("gamma"),
        label = HTML("Posterior probability threshold (\\( \\gamma \\))"),
        min = 80, max = 99, value = 90, step = 1, post = "%"
      ),
      uiOutput(ns("decision_rule"))
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'alpha'", ns("decision_mode")),
      sliderInput(
        ns("alpha"),
        label = HTML("Target Type-I error (\\( \\alpha \\))"),
        value = 10, min = 1, max = 20, step = 1, post = "%"
      ),
      selectInput(
        ns("calibrate_on"),
        "Calibrate Type-I on:",
        choices = c(
          "Point estimate Pr(reject | H₀)" = "point",
          "Upper 95% MC CI (conservative)" = "upper",
          "Lower 95% MC CI (liberal)"      = "lower"
        ),
        selected = "upper"
      ),
      helpText("γ will be calibrated so that the chosen Type-I metric ≈ α (within tolerance) at the least-favourable null.")
    )
  )
}
