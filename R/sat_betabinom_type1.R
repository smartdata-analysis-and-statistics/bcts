#' Estimate Bayesian Type-I Error for a Single-Arm Binomial Trial
#'
#' @description Computes the Bayesian Type-I error for a single-arm binomial trial
#' using a conjugate Beta prior. You can choose between a simulation-based approach
#' or an exact summation over all possible outcomes under the null hypothesis.
#'
#' By default, the null hypothesis assumes \code{p_null = M} (i.e., testing at the boundary),
#' but for a more conservative frequentist setting, you may specify \code{p_null < M}.
#'
#' @param B Integer. Number of simulations (only used if \code{method = "simulate"}).
#' @param n_t Integer. Sample size of the treatment arm.
#' @param M Numeric. Decision threshold (also used as null value if \code{p_null} not specified).
#' @param threshold Numeric. Posterior probability threshold for declaring success, e.g., \code{0.95}.
#' @param prior Prior type: one of:
#'   - \code{"beta"} (default) Use custom Beta(\code{a_base}, \code{b_base}) prior..
#'   - \code{"flat"} for a non-informative Beta(1,1) prior;
#'   - \code{"jeffreys"} for the Jeffreys prior Beta(0.5, 0.5);
#' @param a_base Alpha of Beta prior (only used if \code{prior = "beta"}).
#' @param b_base Beta of Beta prior (only used if \code{prior = "beta"}).
#' @param p_null Optional numeric. True response rate under the null hypothesis (default is \code{M}).
#' Only used when \code{method = "exact"}.
#' @param show_progress Logical. Show progress bar during simulation (only relevant for \code{method = "simulate"}).
#' @param method Character. Either \code{"exact"} (default) or \code{"simulate"}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{estimate}}{Estimated Type-I error rate.}
#'   \item{\code{mc_se}}{Monte Carlo standard error (NA for exact).}
#'   \item{\code{rejections}}{Number of rejections (NA for exact).}
#'   \item{\code{B}}{Number of simulations (NA for exact).}
#' }
#'
#' @examples
#' # Simulated Type-I error at boundary (p_null = M)
#' sat_betabinom_type1(B = 10000, n_t = 40, M = 0.65, threshold = 0.9, method = "simulate")
#'
#' # Exact Type-I error at boundary (p_null = M)
#' sat_betabinom_type1(n_t = 40, M = 0.65, threshold = 0.9, method = "exact")
#'
#' # Exact Type-I error with p_null < M (frequentist-style)
#' sat_betabinom_type1(n_t = 40, M = 0.65, threshold = 0.9, method = "exact", p_null = 0.60)
#'
#' @seealso \code{\link{sat_betabinom_power}}, \code{\link{sat_betabinom_type1_exact}}
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
sat_betabinom_type1 <- function(B = 10000, n_t, M, threshold,
                                prior = c("beta", "flat", "jeffreys"),
                                a_base = 1, b_base = 1,
                                 p_null = NULL,
                                 show_progress = TRUE,
                                 method = c("exact", "simulate")) {
  prior <- match.arg(prior)
  method <- match.arg(method)

  # Convert prior to a_base and b_base
  if (prior == "flat") {
    a_base <- 1
    b_base <- 1
  } else if (prior == "jeffreys") {
    a_base <- 0.5
    b_base <- 0.5
  } else if (prior == "beta") {
    # leave a_base and b_base as provided
  } else {
    stop("Invalid prior. Must be one of 'flat', 'jeffreys', or 'beta'.")
  }

  # Default p_null = M if not specified
  if (is.null(p_null)) p_null <- M

  if (method == "simulate") {
    .Call(`_bcts_sat_betabinom_type1`,
          B, n_t, p_null, M, threshold,
          a_base, b_base, show_progress)
  } else {
    .Call(`_bcts_sat_betabinom_type1_exact`,
          n_t, M, threshold,
          a_base, b_base, p_null)
  }
}
