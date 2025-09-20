# inst/app/R/ui_sat_decision_criteria.R

sat_decision_criteria_ui <- function(id = NULL) {
  ns <- NS(id)

  wellPanel(
    h4("Decision Criteria"),

    sliderInput(
      ns("M_sa"),
      label = HTML("Decision Margin (\\( M \\))"),
      min = 0, max = 100, value = 65, step = 1, post = "%"
    ),

    helpText(
      "Define the threshold on the response rate (M) that must be exceeded to declare success.",
      "Higher values make the design more stringent."
    ),

    radioButtons(
      ns("decision_mode_sa"),
      label = "Threshold Specification Method:",
      choices = c(
        "Specify posterior probability threshold (\\( \\gamma \\))" = "gamma",
        "Specify target Type-I error (\\( \\alpha \\))" = "alpha"
      ),
      selected = "gamma"
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'gamma'", ns("decision_mode_sa")),
      sliderInput(
        ns("gamma_sa"),
        label = HTML("Posterior probability threshold (\\( \\gamma \\))"),
        min = 80, max = 99, value = 90, step = 1, post = "%"
      ),
      uiOutput(ns("decision_rule"))
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'alpha'", ns("decision_mode_sa")),
      sliderInput(
        ns("alpha_sa"),
        label = HTML("Target Type-I error (\\( \\alpha \\))"),
        min = 1, max = 20, value = 10, step = 1, post = "%"
      ),
      selectInput(
        ns("calibrate_on_sa"),
        label = "Calibrate Type-I error based on:",
        choices = c(
          "Point estimate Pr(reject | H₀)"      = "point",
          "Upper 95% MC CI (conservative)"      = "upper",
          "Lower 95% MC CI (liberal)"           = "lower"
        ),
        selected = "upper"
      ),
      helpText(
        "The threshold γ will be calibrated such that the chosen Type-I error metric",
        "approximates α at the least-favorable null scenario."
      )
    )
  )
}
