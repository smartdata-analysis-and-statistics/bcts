#' Plot prior vs current-data weight (ESS)
#'
#' Shows the effective sample size (ESS) contributed by the prior on the
#' control arm, decomposed into Baseline Beta and Borrowed history (if any),
#' versus the observed control sample size.
#'
#' Works for both flat and power priors; borrowed component will be zero if
#' no external data or a0 = 0.
#'
#' @param x A `bayesNI` object with $settings$prior, $settings$n_c,
#'   and $settings$prior_args = list(a_base, b_base, a0, n_0, y_0).
#' @param orientation "h" (horizontal) or "v" (vertical).
#' @param include_labels logical; print ESS numbers on bars.
#' @param ... unused.
#' @return A ggplot object.
#' @importFrom rlang .data
#' @export
plot_prior_weight <- function(x, orientation = c("h","v"),
                              include_labels = TRUE, ...) {
  UseMethod("plot_prior_weight")
}

#' @export
plot_prior_weight.bayesNI <- function(x, orientation = c("h","v"),
                                      include_labels = TRUE, ...) {
  orientation <- match.arg(orientation)

  s  <- x$settings
  pa <- s$prior_args %||% list()

  # Baseline Beta pseudo-counts (default 1,1 if missing)
  a_base <- pa$a_base %||% 1
  b_base <- pa$b_base %||% 1
  ess_baseline <- a_base + b_base

  # Borrowed history (power prior); zero if flat/no external
  a0 <- pa$a0 %||% 0
  n0 <- pa$n_0 %||% 0
  ess_borrow <- a0 * n0

  ess_prior   <- ess_baseline + ess_borrow
  ess_current <- s$n_c

  # Data for stacked prior vs current
  df <- tibble::tibble(
    group     = factor(c("Prior","Prior","Current data"),
                       levels = c("Prior","Current data")),
    component = c("Baseline Beta","Borrowed history","Observed n"),
    ESS       = c(ess_baseline, ess_borrow, ess_current)
  ) |>
    dplyr::filter(.data$ESS > 0 | .data$component == "Observed n")

  # Titles adapt to prior type
  prior_label <- if (identical(s$prior, "power")) "Power prior" else "Baseline-only prior"
  subtitle <- if (identical(s$prior, "power")) {
    sprintf("%s: a0 = %.2f, y0 = %s, n0 = %s; Beta(a,b) = Beta(%.1f, %.1f)\nESS(prior) = %.1f (baseline %.1f + borrowed %.1f); current n = %d",
            prior_label, a0, pa$y_0 %||% "-", n0, a_base, b_base,
            ess_prior, ess_baseline, ess_borrow, ess_current)
  } else {
    sprintf("%s: Beta(a,b) = Beta(%.1f, %.1f)\nESS(prior) = %.1f (baseline only); current n = %d",
            prior_label, a_base, b_base, ess_prior, ess_current)
  }

  pos <- ggplot2::position_stack(vjust = 0.5)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$ESS, fill = .data$component)) +
    ggplot2::geom_col(width = 0.6, colour = "black") +
    { if (isTRUE(include_labels))
      ggplot2::geom_text(ggplot2::aes(label = round(.data$ESS, 1)), position = pos, size = 3.6)
      else ggplot2::geom_blank() } +
    ggplot2::geom_hline(yintercept = ess_current, linetype = "dashed", linewidth = 0.4) +
    ggplot2::labs(
      title    = "Relative weight of prior vs current control likelihood",
      subtitle = subtitle,
      x = NULL, y = "Effective sample size (ESS)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (orientation == "h") p <- p + ggplot2::coord_flip()
  p
}

# tiny helper
`%||%` <- function(x, y) if (is.null(x)) y else x
