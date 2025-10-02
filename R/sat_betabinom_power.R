#' Estimate Bayesian power for a single-arm binomial trial
#'
#' @description Computes the Bayesian power for a single-arm binomial trial
#' using either simulation or exact analytic summation over all possible outcomes.
#'
#' @param B Number of trial simulations (used only if \code{method = "simulate"}).
#' @param p_t True response probability.
#' @param n_t Sample size.
#' @param M Threshold for success (e.g. 0.6).
#' @param threshold Posterior probability threshold (e.g. 0.95).
#' @param prior Prior type: one of:
#'   - \code{"flat"} for a non-informative Beta(1,1) prior;
#'   - \code{"jeffreys"} for the Jeffreys prior Beta(0.5, 0.5);
#'   - \code{"beta"} for a custom Beta(\code{a_base}, \code{b_base}) prior.
#' @param a_base Alpha of Beta prior (only used if \code{prior = "beta"}).
#' @param b_base Beta of Beta prior (only used if \code{prior = "beta"}).
#' @param show_progress Show progress bar? (only relevant for simulation).
#' @param method Either \code{"simulate"} (default) or \code{"exact"}.
#'
#' @return A list with \code{estimate} (power), \code{mc_se}, \code{successes}, and \code{B}.
#'
#' @examples
#' sat_betabinom_power(
#'   B = 1000, p_t = 0.75, n_t = 35, M = 0.60,
#'   threshold = 0.95, prior = "flat", method = "simulate"
#' )
#'
#' sat_betabinom_power(
#'   p_t = 0.75, n_t = 35, M = 0.60,
#'   threshold = 0.95, prior = "flat", method = "exact"
#' )
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @export
sat_betabinom_power <- function(B = 10000, p_t, n_t, M, threshold,
                                prior = c("flat", "jeffreys", "beta"),
                                a_base = 1, b_base = 1,
                                 show_progress = TRUE,
                                 method = c("exact", "simulate")) {
  prior <- match.arg(prior)
  method <- match.arg(method)

  if (method == "simulate") {
    .Call(`_bcts_sat_betabinom_power`,
          B, p_t, n_t, M, threshold,
          prior, a_base, b_base, show_progress)
  } else {
    .Call(`_bcts_sat_betabinom_power_exact`,
          p_t, n_t, M, threshold,
          prior, a_base, b_base)
  }
}
