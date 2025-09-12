#' Run a full Bayesian NI analysis (Beta–Binomial, conjugate)
#'
#' Calibrates the posterior cutoff \eqn{\gamma} to control Type-I error at a target
#' \code{alpha}, and then estimates the corresponding Bayesian power for given trial
#' parameters. Uses conjugate Beta–Binomial updates (no MCMC).
#'
#' @param alpha Target Type-I error rate (numeric in `(0,1)`).
#' @param p_c,p_t True response probabilities in control and treatment (numeric in `[0,1]`).
#' @param n_c,n_t Sample sizes in control and treatment arms (integers).
#' @param NI_margin Non-inferiority margin on the risk-difference scale (numeric).
#' @param prior Character. Either `"flat"` (Beta(1,1) both arms) or `"power"`
#'   (power prior on control); see `prior_args`.
#' @param prior_args List of hyperparameters when `prior = "power"`:
#'   - `a0`: discount factor in `[0,1]`
#'   - `y_0`, `n_0`: historical control responders and sample size
#'   - `a_base`, `b_base`: baseline Beta parameters (defaults `1`, `1`)
#' @param B_cal Number of simulated trials for calibration of `gamma`.
#' @param B_power Number of simulated trials for power estimation.
#' @param n_draws Posterior Monte Carlo draws per trial for evaluating `Pr(NI)`.
#' @param tol Absolute tolerance for calibration convergence.
#' @param seed RNG seed for reproducibility.
#' @param show_progress Logical. Show text progress bars during simulation.
#' @param verbose Logical. Print iteration messages during calibration.
#'
#' @return An object of class `"bayesNI"` with components:
#' - `calibration`: list from [bayesNI_calibrate_betaBinom_conj()]
#' - `power`: estimated Bayesian power at calibrated `gamma`
#' - `settings`: list of all input parameters (echoed)
#'
#' @seealso [bayesNI_calibrate_betaBinom_conj()],
#'   [bayesNI_power_betaBinom_conj()],
#'   [bayesNI_type1_betaBinom_conj()]
#' @family conjugate-NI
#'
#' @examples
#' \donttest{
#' res <- bayesNI_run_betaBinom_conj(
#'   alpha = 0.10, p_c = .85, p_t = .85,
#'   n_c = 29, n_t = 29, NI_margin = -0.20,
#'   prior = "flat", B_cal = 1000, B_power = 500
#' )
#' res$calibration$gamma
#' res$power
#' }
#'
#' @export
bayesNI_run_betaBinom_conj <- function(alpha = 0.10,
                              p_c = 0.85, p_t = 0.85,
                              n_c = 29, n_t = 29,
                              NI_margin = -0.20,
                              prior = c("flat","power"),
                              prior_args = list(),
                              B_cal = 2000, B_power = 1000,
                              n_draws = 2000,
                              tol = 0.002,
                              seed = 123,
                              show_progress = FALSE, verbose = TRUE) {
  prior <- match.arg(prior)

  # 1) Calibrate gamma using conjugate machinery
  cal <- bayesNI_calibrate_betaBinom_conj(alpha = alpha, p_c = p_c, M = NI_margin,
                              n_c = n_c, n_t = n_t,
                              prior = prior, prior_args = prior_args,
                              B_cal = B_cal, n_draws = n_draws,
                              tol = tol, seed = seed,
                              show_progress = show_progress, verbose = verbose)

  # 2) Simulate Bayesian "power" at calibrated gamma
  power <- bayesNI_power_betaBinom_conj(B = B_power, p_c = p_c, p_t = p_t,
                                     n_c = n_c, n_t = n_t, M = NI_margin,
                                     threshold = cal$gamma,
                                     prior = prior, prior_args = prior_args,
                                     n_draws = n_draws,
                                     seed = seed, show_progress = show_progress)

  out <- list(
    calibration = cal,
    power = power,
    settings = list(
      alpha = alpha, p_c = p_c, p_t = p_t,
      n_c = n_c, n_t = n_t, NI_margin = NI_margin,
      prior = prior, prior_args = prior_args,
      B_cal = B_cal, B_power = B_power, tol = tol, seed = seed,
      engine = "conjugate", n_draws = n_draws
    )
  )
  class(out) <- "bayesNI"
  out
}
