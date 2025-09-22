# modules/mod_narrative_ui
message("Loaded mod_narrative_ui.R")

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
          "This design simulates a two-arm trial to %s. The threshold γ is calibrated to control the Type-I error at %s.",
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
        "flat" = "a non-informative flat prior for the response rate in each arm",
        "power" = {
          a0   <- prior_args()$a0   %||% 0
          n_0  <- prior_args()$n_0  %||% 1
          y_0  <- prior_args()$y_0  %||% 0
          eff_n0 <- a0 * n_0

          sprintf(
            "a power prior incorporating %s/%s historical responders. The historical data is down-weighted using a factor of %.0f%%, resulting in an effective contribution equivalent to %.1f patients",
            fmt_int(y_0),
            fmt_int(n_0),
            100 * a0,
            eff_n0
          )
        }
      )

      # 6. Final narrative
      HTML(sprintf(
        "<p>%s</p>
         <p>%s</p>
         <p>The analysis uses %s.</p>
         <p>Posterior probabilities are computed using %s draws per simulated trial, across %s simulated trials.</p>",
        design_txt,
        arm_txt,
        prior_type,
        fmt_int(ndraws()), fmt_int(B())
      ))
    })
  })
}
