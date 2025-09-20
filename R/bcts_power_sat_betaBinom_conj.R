#' Simulate Bayesian power for a single-arm binomial trial (via Rcpp)
#'
#' @param B Number of trial simulations.
#' @param p_t True response probability.
#' @param n_t Sample size.
#' @param M Threshold for success (e.g. 0.6).
#' @param threshold Posterior probability threshold (e.g. 0.95).
#' @param prior Prior type: "flat" or "beta".
#' @param a_base Alpha of Beta prior (only used if prior = "beta").
#' @param b_base Beta of Beta prior (only used if prior = "beta").
#' @param n_draws Number of posterior samples per trial.
#' @param show_progress Show progress bar?
#'
#' @return A list with `estimate` (power), `mc_se`, `successes`, and `B`.
#'
#' @export
singlearm_beta_power <- function(B, p_t, n_t, M, threshold,
                                 prior = "flat", a_base = 1, b_base = 1,
                                 n_draws = 2000, show_progress = TRUE) {
  .Call(`_bcts_singlearm_beta_power`, B, p_t, n_t, M, threshold,
        prior, a_base, b_base, n_draws, show_progress)
}
