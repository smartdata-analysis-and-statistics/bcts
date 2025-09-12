#' Bayesian NI power via conjugate Beta–Binomial
#'
#' Monte Carlo estimate of Bayesian non-inferiority (NI) \emph{power} using
#' conjugate Beta–Binomial updates (no MCMC). Each replicate simulates a trial,
#' computes the posterior probability of NI, and compares it to a fixed
#' threshold \code{threshold}.
#'
#' @param B Integer. Number of simulated trials.
#' @param p_c,p_t Numeric in \[0,1\]. True response probabilities in control and treatment.
#' @param n_c,n_t Integers. Sample sizes in control and treatment.
#' @param M Numeric. NI margin on the risk-difference scale (e.g., \code{-0.20}).
#' @param threshold Numeric in (0,1). Posterior probability cutoff (gamma) for declaring NI.
#' @param prior Either `"flat"` (Beta(1,1) both arms) or `"power"` (power prior on control).
#' @param prior_args List of hyperparameters when `prior = "power"`:
#'   - `a0`: discount factor in `[0,1]`
#'   - `y_0`, `n_0`: historical control responders and sample size
#'   - `a_base`, `b_base`: baseline Beta parameters (defaults `1`, `1`)
#' @param n_draws Integer. Posterior Monte Carlo draws per trial (for \code{Pr(NI)} evaluation).
#' @param seed Optional integer RNG seed (reproducibility of the outer loop).
#' @param show_progress Logical. Show a text progress bar.
#'
#' @return Numeric scalar: estimated Bayesian power (i.e., `Pr(declare NI)`).
#'
#' @examples
#' \donttest{
#' bayesNI_power_betaBinom_conj(
#'   B = 200, p_c = .85, p_t = .85, n_c = 29, n_t = 29,
#'   M = -0.20, threshold = 0.91, prior = "flat", n_draws = 1000
#' )
#' }
#'
#' @seealso \code{\link{bayesNI_calibrate_betaBinom_conj}}, \code{\link{bayesNI_type1_betaBinom_conj}}
#' @family conjugate-NI
#' @export
bayesNI_power_betaBinom_conj <- function(B = 1000, p_c, p_t, n_c, n_t, M,
                                      threshold, prior = c("flat","power"),
                                      prior_args = list(),
                                      n_draws = 2000, seed = NULL,
                                      show_progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  decide <- logical(B)
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  for (b in seq_len(B)) {
    out <- bayesNI_trial_betaBinom_conj(p_c = p_c, p_t = p_t, n_c = n_c, n_t = n_t, M = M,
                               prior = prior, prior_args = prior_args,
                               n_draws = n_draws)
    decide[b] <- (out$summary["post_prob_NI"] >= threshold)
    if (show_progress) utils::setTxtProgressBar(pb, b)
  }
  mean(decide)
}

#' Empirical Type-I error for Bayesian NI (conjugate)
#'
#' Estimates the Bayesian NI \emph{Type-I error} at a fixed posterior threshold
#' \code{threshold} under the least-favourable null (LFN) \code{p_t = p_c + M}.
#' Uses conjugate Beta–Binomial updates (no MCMC).
#'
#' @param B Integer. Number of simulated trials.
#' @param p_c Numeric in \[0,1\]. True control response probability.
#' @param M Numeric. NI margin on the risk-difference scale.
#' @param n_c,n_t Integers. Sample sizes in control and treatment.
#' @param threshold Numeric in (0,1). Posterior probability cutoff (gamma).
#' @param prior Character. \code{"flat"} or \code{"power"}; see Details and \code{prior_args}.
#' @param prior_args List for power prior on control:
#'   - `a0`, `y_0`, `n_0`, `a_base`, `b_base`
#' @param n_draws Posterior Monte Carlo draws per trial for evaluating `Pr(NI)`.
#' @param seeds Optional integer vector of length `B` (common random numbers).
#'   If `NULL`, a default set is generated.
#' @param show_progress Show a text progress bar.
#'
#' @return A list with:
#' - `type1`: estimated Type-I error
#' - `mc_se`: Monte Carlo standard error
#' - `mc_ci`: approximate 95% MC confidence interval
#' - `seeds`: vector of seeds actually used
#' - `settings`: echo of key inputs, including the LFN `p_t`
#'
#' @examples
#' \donttest{
#' bayesNI_type1_betaBinom_conj(
#'   B = 500, p_c = .85, M = -0.20, n_c = 29, n_t = 29,
#'   threshold = 0.91, prior = "flat", n_draws = 1000
#' )
#' }
#'
#' @seealso [bayesNI_calibrate_betaBinom_conj()], [bayesNI_power_betaBinom_conj()]
#' @family conjugate-NI
#' @export
bayesNI_type1_betaBinom_conj <- function(B = 2000, p_c, M, n_c, n_t,
                                      threshold, prior = c("flat","power"),
                                      prior_args = list(),
                                      n_draws = 2000,
                                      seeds = NULL, show_progress = TRUE) {
  # LFC: p_t = p_c + M
  p_t <- p_c + M
  if (p_t < 0 || p_t > 1) stop("Least-favourable null p_t = p_c + M outside [0,1].")

  if (is.null(seeds)) { set.seed(1L); seeds <- sample.int(1e9, B) }
  decide <- logical(B)

  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  for (b in seq_len(B)) {
    set.seed(seeds[b])
    out <- bayesNI_trial_betaBinom_conj(p_c = p_c, p_t = p_t, n_c = n_c, n_t = n_t, M = M,
                               prior = prior, prior_args = prior_args,
                               n_draws = n_draws)
    decide[b] <- (out$summary["post_prob_NI"] >= threshold)
    if (show_progress) utils::setTxtProgressBar(pb, b)
  }

  type1 <- mean(decide)
  se    <- sqrt(type1 * (1 - type1) / B)
  ci    <- type1 + c(-1, 1) * 1.96 * se

  list(type1 = type1, mc_se = se, mc_ci = ci, seeds = seeds,
       settings = list(B = B, p_c = p_c, p_t = p_t, n_c = n_c, n_t = n_t,
                       M = M, threshold = threshold, prior = match.arg(prior)))
}

#' Calibrate the posterior threshold \eqn{\gamma} (conjugate)
#'
#' Finds the posterior probability cutoff \eqn{\gamma} such that the Bayesian NI
#' decision rule \code{Pr(theta_t - theta_c > M | data) >= gamma} achieves a
#' target Type-I error \code{alpha} under the least-favourable null, using
#' conjugate Beta–Binomial updates (no MCMC). A bisection search is performed
#' over \code{[lower, upper]} with caching and common random numbers.
#'
#' @param alpha Target Type-I error in `(0,1)`, e.g., `0.10`.
#' @param p_c Control response probability under the null (in `[0,1]`).
#' @param M NI margin on the risk-difference scale.
#' @param n_c,n_t Sample sizes in control and treatment (integers).
#' @param prior Either `"flat"` or `"power"`. See `prior_args`.
#' @param prior_args List for power prior on control:
#'   - `a0`, `y_0`, `n_0`, `a_base`, `b_base`
#' @param B_cal Number of simulated trials per Type-I evaluation.
#' @param lower,upper Initial bracketing interval for `gamma`, values in `(0,1)`.
#' @param n_draws Posterior Monte Carlo draws per trial.
#' @param tol Absolute tolerance for `|TypeI - alpha|`. If `NULL`, a Monte
#'   Carlo–driven default is used.
#' @param maxit Maximum number of bisection iterations.
#' @param seed RNG seed for generating common random numbers across thresholds.
#' @param show_progress Show a progress bar inside Type-I evaluations.
#' @param verbose Print bracketing/iteration messages.
#' @param digits Rounding used in cache keys (affects only speed).
#'
#' @return A list with:
#' - `gamma`: calibrated posterior threshold
#' - `type1`: estimated Type-I error at `gamma`
#' - `bracket`: final bracketing interval
#' - `iters`: number of bisection iterations used
#' - `B_cal`, `alpha`: echoed inputs
#' - `trace`: data frame with the search history
#'
#' @examples
#' \donttest{
#' bayesNI_calibrate_betaBinom_conj(
#'   alpha = 0.10, p_c = .85, M = -0.20, n_c = 29, n_t = 29,
#'   prior = "flat", B_cal = 1000, n_draws = 1000, seed = 11
#' )
#' }
#'
#' @seealso [bayesNI_type1_betaBinom_conj()], [bayesNI_power_betaBinom_conj()]
#' @family conjugate-NI
#' @export
bayesNI_calibrate_betaBinom_conj <- function(alpha = 0.10, p_c, M, n_c, n_t,
                                 prior = c("flat","power"),
                                 prior_args = list(),
                                 B_cal = 2000, lower = 0.80, upper = 0.999,
                                 n_draws = 2000,
                                 tol = 0.002, maxit = 30,
                                 seed = 11L, show_progress = FALSE,
                                 verbose = TRUE, digits = 4) {

  prior <- match.arg(prior)
  set.seed(seed)
  seeds <- sample.int(1e9, B_cal)

  if (is.null(tol)) {
    tol <- 1.25 * sqrt(alpha * (1 - alpha) / B_cal)
    if (verbose) message(sprintf("Auto tol set to %.4f based on B_cal=%d", tol, B_cal))
  }

  cache <- new.env(parent = emptyenv())
  keyfun <- function(thr) paste0("thr_", format(round(thr, digits), nsmall = digits))

  type1_at <- function(thr) {
    key <- keyfun(thr)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    est <- bayesNI_type1_betaBinom_conj(B = B_cal, p_c = p_c, M = M, n_c = n_c, n_t = n_t,
                                     threshold = thr, prior = prior, prior_args = prior_args,
                                     n_draws = n_draws, seeds = seeds,
                                     show_progress = show_progress)$type1
    assign(key, est, envir = cache)
    est
  }

  t_lo <- type1_at(lower)
  t_hi <- type1_at(upper)
  if (verbose) {
    message(sprintf("Init: lower=%.3f -> type1=%.4f | upper=%.3f -> type1=%.4f | target alpha=%.4f",
                    lower, t_lo, upper, t_hi, alpha))
  }
  if (!(t_lo >= alpha && t_hi <= alpha)) {
    warning(sprintf("Alpha not bracketed: type1(lower)=%.4f, type1(upper)=%.4f", t_lo, t_hi))
  }

  lo <- lower; hi <- upper
  trace <- data.frame(iter = integer(0),
                      gamma_try = numeric(0),
                      type1 = numeric(0),
                      lo = numeric(0),
                      hi = numeric(0),
                      diff = numeric(0))

  for (i in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    t_mid <- type1_at(mid)
    diff  <- t_mid - alpha
    trace <- rbind(trace, data.frame(iter = i, gamma_try = mid, type1 = t_mid,
                                     lo = lo, hi = hi, diff = diff))
    if (verbose) message(sprintf("Iter %02d: gamma=%.4f -> type1=%.4f (diff=%.4f)", i, mid, t_mid, diff))
    if (abs(diff) < tol) {
      return(list(gamma = mid, type1 = t_mid,
                  bracket = c(lower = lo, upper = hi),
                  iters = i, B_cal = B_cal, alpha = alpha,
                  trace = trace))
    }
    if (t_mid > alpha) lo <- mid else hi <- mid
  }

  mid <- (lo + hi) / 2
  t_mid <- type1_at(mid)
  trace <- rbind(trace, data.frame(iter = maxit + 1, gamma_try = mid, type1 = t_mid,
                                   lo = lo, hi = hi, diff = t_mid - alpha))
  list(gamma = mid, type1 = t_mid,
       bracket = c(lower = lo, upper = hi),
       iters = maxit, B_cal = B_cal, alpha = alpha, trace = trace)
}
