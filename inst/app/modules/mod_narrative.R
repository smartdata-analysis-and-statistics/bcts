
# ---- UI ----
mod_narrative_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Narrative summary"),
    htmlOutput(ns("narrative_text"))
  )
}

# ---- Server ----
mod_narrative_server <- function(
    id,
    pt, nt,
    pc, nc,
    M,
    prior, prior_args,
    decision_mode, gamma, alpha,
    B, ndraws, seed
) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(x, y) if (is.null(x)) y else x
    fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
    fmt_int <- function(x) formatC(as.integer(x), big.mark = ",", format = "d")

    output$narrative_text <- renderUI({

      # 1. Margin description
      margin_txt <- if (M() > 0) {
        sprintf("demonstrate that the active treatment is superior to control (Δ = +%s)", fmt_pct(M()))
      } else {
        sprintf("demonstrate non-inferiority of the active treatment, allowing a margin of Δ = %s", fmt_pct(M()))
      }

      # 2. Trial design description
      design_txt <- if (decision_mode() == "gamma") {
        sprintf(
          "This design simulates a two-arm trial to %s. The trial is declared successful if the posterior probability that the true treatment effect exceeds Δ is at least γ = %s.",
          margin_txt, fmt_pct(gamma())
        )
      } else {
        sprintf(
          "This design simulates a two-arm trial to %s, using a one-sided significance level of %s.",
          margin_txt, fmt_pct(alpha())
        )
      }

      # 3. Treatment and control arm setup
      arm_txt <- sprintf(
        "The treatment arm includes %s patients with an assumed true response rate of %.1f%%.
         The control arm includes %s patients with an assumed true response rate of %.1f%%.",
        fmt_int(nt()), 100 * pt(),
        fmt_int(nc()), 100 * pc()
      )

      # 4. Prior description
      prior_type <- switch(
        prior(),
        "flat" = "a non-informative flat prior (Beta(1,1))",
        "power" = sprintf(
          "a power prior (discount a₀ = %.0f%%) incorporating %s/%s historical responders and a baseline Beta prior (%.1f, %.1f)",
          100 * prior_args()$a0,
          fmt_int(prior_args()$y_0),
          fmt_int(prior_args()$n_0),
          prior_args()$a_base %||% 1,
          prior_args()$b_base %||% 1
        )
      )

      # 5. Decision logic
      decision_txt <- if (decision_mode() == "gamma") {
        sprintf("The trial is declared successful if the posterior probability that the true treatment effect exceeds Δ is at least γ = %s.", fmt_pct(gamma()))
      } else {
        sprintf("γ is calibrated to achieve a target Type-I error of %s.", fmt_pct(alpha()))
      }

      # 6. Final narrative
      HTML(sprintf(
        "<p>%s</p>
         <p>%s</p>
         <p>%s</p>
         <p>The analysis uses %s.</p>
         <p>Posterior probabilities are computed using %s draws per simulated trial, across %s simulated trials.</p>
         <p>Seed = %s (for reproducibility).</p>",
        design_txt,
        arm_txt,
        decision_txt,
        prior_type,
        fmt_int(ndraws()), fmt_int(B()),
        fmt_int(seed())
      ))
    })
  })
}
