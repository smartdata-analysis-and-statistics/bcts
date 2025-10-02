#' Calibrate posterior cutoff gamma and estimate power (Single-Arm, Beta–Binomial)
#'
#' Calibrates the posterior threshold `gamma` to control Type-I error at level `alpha`,
#' and estimates Bayesian power using a conjugate Beta–Binomial model.
#'
#' @param alpha Target Type-I error rate (numeric in `(0,1)`).
#' @param p_t True response probability in treatment group (numeric in `[0,1]`).
#' @param n_t Sample size in treatment group (integer).
#' @param M Decision threshold for theta (numeric in `[0,1]`).
#' @param a_base,b_base Numeric. Hyperparameters of the Beta prior (used when `prior = "beta"`).
#' @param tol Absolute tolerance for calibration convergence.
#' @param seed RNG seed.
#' @param show_progress Logical. Whether to display progress bars.
#'
#' @return A list with components:
#'   - `gamma`: Calibrated posterior probability threshold.
#'   - `type1`: Estimated Type-I error at `gamma`.
#'   - `power`: Estimated Bayesian power at `gamma`.
#'   - `settings`: Echo of input parameters.
#'
#' @examples
#' \donttest{
#' res <- sat_beta_calibrate(
#'   alpha = 0.10, p_t = 0.80, n_t = 40, M = 0.65,
#'   prior = "flat"
#' )
#' res$gamma
#' res$power
#' }
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
sat_betabinom_calibrate <- function(alpha = 0.10,
                              p_t = 0.80,
                              n_t = 40,
                              M = 0.65,
                              prior = c("flat", "jeffreys", "beta"),
                              a_base = 1, b_base = 1,
                              tol = 0.002,
                              seed = 123,
                              show_progress = FALSE) {
  prior <- match.arg(prior)

  cal <- sat_betabinom_find_gamma(
    alpha = alpha, p_c = M, M = M,
    n_t = n_t,
    prior = prior,
    prior_args = list(a_base = a_base, b_base = b_base),
    seed = seed,
    show_progress = show_progress,
    progress_fun = pf,
    calibrate_on = calibrate_on
  )


  power <- sat_betabinom_power(
    p_t = p_t, n_t = n_t, M = M,
    threshold = cal$gamma,
    prior = prior, a_base = a_base, b_base = b_base,
    show_progress = show_progress,
    method = "exact"
  )

  type1 <- sat_betabinom_type1(
    n_t = n_t, M = M,
    threshold = cal$gamma,
    prior = prior, a_base = a_base, b_base = b_base,
    show_progress = show_progress,
    method = "exact"
  )

  out <- list(
    gamma = cal$gamma,
    type1 = type1,
    power = power,
    settings = list(
      alpha = alpha, p_t = p_t, n_t = n_t, M = M,
      prior = prior, a_base = a_base, b_base = b_base,
      B_cal = B_cal, B_power = B_power,
      tol = tol, seed = seed
    )
  )
  class(out) <- "bcts_sat"
  return(out)
}

#' Calibrate posterior threshold gamma to control Type-I error (exact, SAT)
#'
#' For single-arm trials with Beta–Binomial conjugate updates, computes the Type-I error
#' exactly for all possible posterior thresholds, and selects the largest gamma such that
#' Type-I error is below or equal to a target alpha.
#'
#' @param alpha Target Type-I error rate (e.g. 0.10)
#' @param n_t Sample size in treatment arm
#' @param M Decision margin (null threshold)
#' @param a_base,b_base Prior Beta(a,b) parameters (default: 1,1)
#'
#' @return An object of class `"bcts_calibration"` with components: `gamma`, `type1`, `trace`, and `settings`.
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
sat_betabinom_calibrate_exact <- function(alpha = 0.10,
                                          n_t,
                                          M,
                                          a_base = 1,
                                          b_base = 1,
                                          grid_n = 1000) {

  # Step 1: Grid of gamma values (search over [0,1])
  gamma_vals <- seq(1, 0, length.out = grid_n)  # from most to least conservative

  # Step 2: Evaluate type-I error at each gamma
  type1_vals <- vapply(gamma_vals, function(g) {
    .Call(`_bcts_sat_betabinom_type1_exact`,
          n_t, M, g, a_base, b_base, -1.0)
  }, numeric(1))

  # Step 3: Build trace table
  trace <- data.frame(
    gamma = gamma_vals,
    type1 = type1_vals
  )

  # Step 4: Select gamma where type-I is just below alpha
  eligible <- trace$type1 <= alpha
  if (!any(eligible)) {
    stop("Calibration failed: Type-I error is above target alpha at all gamma values.")
  }

  best_idx <- which.max(trace$gamma[eligible])  # most liberal gamma ≤ alpha
  best_gamma <- trace$gamma[eligible][best_idx]
  best_type1 <- trace$type1[eligible][best_idx]

  out <- list(
    gamma = best_gamma,
    type1 = best_type1,
    trace = trace
  )
  class(out) <- "bcts_calibration"

  return(out)
}
