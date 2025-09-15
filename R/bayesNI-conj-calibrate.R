#' Bayesian trial power estimation with Beta–Binomial conjugate model
#'
#' Uses Monte Carlo simulation with conjugate Beta–Binomial updates
#' (no MCMC) to estimate the probability of declaring trial success
#' under user-specified assumptions. A trial is declared successful
#' if \eqn{Pr(\theta_t - \theta_c > M | data) \ge \gamma}.
#'
#' @param B Integer. Number of simulated trials.
#' @param p_c,p_t Numeric in \[0,1\]. True response probabilities in control and treatment.
#' @param n_c,n_t Integers. Sample sizes in control and treatment.
#' @param M Numeric. Margin on the risk-difference scale.
#'   - Use negative values for non-inferiority,
#'   - Zero for equivalence,
#'   - Positive for superiority.
#' @param threshold Numeric in (0,1). Posterior probability cutoff \eqn{\gamma}.
#' @param prior Either `"flat"` (Beta(1,1) both arms) or `"power"` (power prior on control).
#' @param prior_args List of hyperparameters when `prior = "power"`:
#'   - `a0`: discount factor in `[0,1]`
#'   - `y_0`, `n_0`: historical control responders and sample size
#'   - `a_base`, `b_base`: baseline Beta parameters (defaults `1`, `1`)
#' @param n_draws Integer. Posterior Monte Carlo draws per trial (for \eqn{Pr()} evaluation).
#' @param seed Optional integer RNG seed.
#' @param show_progress Logical. Show a text progress bar.
#' @param conf.level Numeric in (0,1). Confidence level for binomial CI.
#' @param ci_method Character. CI method passed to [binom::binom.confint()].
#'
#' @return A list with:
#'   - `estimate`: estimated success probability (power)
#'   - `ci`: matrix with lower/upper bounds for the CI
#'   - `B`: number of simulated trials
#'   - `decisions`: logical vector of trial-level NI declarations
#'
#' @examples
#' res <- bcts_power_betaBinom_conj(
#'   B = 100, p_c = 0.85, p_t = 0.85, n_c = 29, n_t = 29,
#'   M = -0.20, threshold = 0.90, prior = "flat", n_draws = 200, seed = 123
#' )
#' res$estimate
#' res$ci
#'
#' @seealso [bcts_type1_betaBinom_conj()], [bcts_calibrate_betaBinom_conj()]
#' @export
bcts_power_betaBinom_conj <- function(B = 1000, p_c, p_t, n_c, n_t, M,
                                         threshold,
                                         prior = c("flat","power"),
                                         prior_args = list(),
                                         n_draws = 2000, seed = NULL,
                                         show_progress = TRUE,
                                         conf.level = 0.95,
                                         ci_method = "wilson") {
  prior <- match.arg(prior)
  if (!is.null(seed)) set.seed(seed)

  decide <- logical(B)
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  for (b in seq_len(B)) {
    out <- bayesNI_trial_betaBinom_conj(
      p_c = p_c, p_t = p_t, n_c = n_c, n_t = n_t, M = M,
      prior = prior, prior_args = prior_args,
      n_draws = n_draws
    )
    decide[b] <- (out$summary["post_prob_NI"] >= threshold)
    if (show_progress) utils::setTxtProgressBar(pb, b)
  }
  k <- sum(decide)
  p <- k / B
  ci <- binom::binom.confint(x=k, n=B, conf.level = conf.level, methods = ci_method)

  list(estimate = p,
       ci_lower = ci["lower"],
       ci_upper = ci["upper"],
       B = B, successes = k)
}

#' Empirical Type-I Error for Bayesian Trial (Beta–Binomial, conjugate)
#'
#' Estimates the Type-I error rate for a Bayesian non-inferiority (or
#' superiority) test using conjugate Beta–Binomial updates. The test declares
#' success if the posterior probability \eqn{Pr(\theta_t - \theta_c > M | data)}
#' exceeds a fixed threshold `threshold`. Simulation is conducted under the
#' least-favourable null (LFN), where \eqn{p_t = p_c + M}.
#'
#' @param B Integer. Number of simulated trials.
#' @param p_c Numeric in \[0,1\]. True response probability in the control arm.
#' @param M Numeric. Margin on the risk-difference scale. The LFN is defined as
#'   \eqn{p_t = p_c + M}.
#' @param n_c,n_t Integers. Sample sizes in control and treatment arms.
#' @param threshold Numeric in (0,1). Posterior probability cutoff (gamma) for
#'   declaring success.
#' @param prior Character. Either `"flat"` (Beta(1,1) prior for both arms) or
#'   `"power"` (power prior on control arm); see `prior_args`.
#' @param prior_args List of hyperparameters if `prior = "power"`, including:
#'   - `a0`: discount factor in \[0,1\]
#'   - `y_0`, `n_0`: historical control responders and sample size
#'   - `a_base`, `b_base`: baseline Beta prior parameters
#' @param n_draws Integer. Posterior Monte Carlo draws per simulated trial.
#' @param seeds Optional integer vector of length `B`. If `NULL`, a default set
#'   of random seeds is generated.
#' @param show_progress Logical. Display a text progress bar during simulation.
#' @param ci_method Character. Method for binomial confidence interval;
#'   `"wilson"` (default) or `"exact"`. Requires the \pkg{binom} package if
#'   `"exact"` is used.
#' @param conf.level Numeric in (0,1). Confidence level for the binomial CI
#'   (default 0.95).
#'
#' @return A list with:
#' - `type1`: estimated Type-I error
#' - `mc_se`: Monte Carlo standard error
#' - `mc_ci`: approximate 95% MC confidence interval
#' - `seeds`: vector of seeds actually used
#' - `settings`: echo of key inputs, including the LFN `p_t`
#'
#' @examples
#' # Quick example (small B and n_draws for speed in testing)
#' bcts_type1_betaBinom_conj(
#'   B         = 50,
#'   p_c       = 0.85,
#'   M         = -0.20,
#'   n_c       = 29,
#'   n_t       = 29,
#'   threshold = 0.90,
#'   prior     = "flat",
#'   n_draws   = 200
#' )
#'
#' @seealso [bcts_calibrate_betaBinom_conj()], [bcts_power_betaBinom_conj()]
#' @family conjugate-BetaBinom
#' @export
bcts_type1_betaBinom_conj <- function(B = 2000, p_c, M, n_c, n_t,
                                      threshold, prior = c("flat","power"),
                                      prior_args = list(),
                                      n_draws = 2000,
                                      seeds = NULL,
                                      show_progress = TRUE,
                                      ci_method = c("wilson","exact"),
                                      conf.level = 0.95) {
  # match/validate args
  prior     <- match.arg(prior)
  ci_method <- match.arg(ci_method)

  # LFN: p_t = p_c + M
  p_t <- p_c + M
  if (p_t < 0 || p_t > 1) stop("Least-favourable null p_t = p_c + M outside [0,1].")

  # seeds
  if (is.null(seeds)) { set.seed(1L); seeds <- sample.int(1e9, B) }

  decide <- logical(B)

  # progress bar
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  for (b in seq_len(B)) {
    set.seed(seeds[b])
    out <- bayesNI_trial_betaBinom_conj(
      p_c = p_c, p_t = p_t,
      n_c = n_c, n_t = n_t, M = M,
      prior = prior, prior_args = prior_args,
      n_draws = n_draws
    )
    decide[b] <- unname(out$summary["post_prob_NI"]) >= threshold
    if (show_progress) utils::setTxtProgressBar(pb, b)
  }

  # summary + CI
  k <- sum(decide)
  p <- k / B
  ci <- binom::binom.confint(x = k, n = B,
                             conf.level = conf.level,
                             methods = ci_method)
  mc_se <- sqrt(p * (1 - p) / B)

  list(
    estimate   = p,
    ci_lower   = ci$lower,
    ci_upper   = ci$upper,
    mc_se      = mc_se,
    B          = B,
    successes  = k,
    seeds      = seeds,
    settings   = list(
      B = B, p_c = p_c, p_t = p_t, n_c = n_c, n_t = n_t,
      M = M, threshold = threshold, prior = prior,
      n_draws = n_draws, conf.level = conf.level, ci_method = ci_method
    )
  )
}

#' Calibrate the posterior threshold gamma (Beta–Binomial, conjugate)
#'
#' Finds the posterior probability cutoff \eqn{\gamma} such that the Bayesian
#' decision rule \code{Pr(theta_t - theta_c > M | data) >= gamma} achieves a
#' target Type-I error \code{alpha} under the least-favourable null, using
#' conjugate Beta–Binomial updates (no MCMC). A bisection search is performed
#' over \code{[lower, upper]} with caching and common random numbers.
#'
#' @param alpha Target Type-I error in `(0,1)`, e.g., `0.10`.
#' @param p_c Control-arm response probability under the null (in `[0,1]`).
#' @param M Margin on the risk-difference scale (negative for NI; nonnegative for superiority).
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
#' @param progress_fun Optional function called once per bisection iteration:
#'   `progress_fun(iter, maxit)`. In Shiny, pass a wrapper that calls
#'   `incProgress(1/maxit, detail = sprintf("Iteration %d of %d", iter, maxit))`.
#' @param calibrate_on Character string indicating which Type-I error quantity
#'   to match during calibration. One of:
#'   * `"point"` — use the Monte Carlo **point estimate** of the Type-I error
#'   * `"upper"` — use the **upper bound** of the 95% Monte Carlo CI
#'   * `"lower"` — use the **lower bound** of the 95% Monte Carlo CI
#'   Default is `"point"`.
#'
#' @return A list with:
#' - `gamma`: calibrated posterior threshold
#' - `type1`: estimated Type-I error at `gamma`
#' - `mc_se`: Monte Carlo standard error at `gamma`
#' - `mc_ci`: approximate 95% MC confidence interval at `gamma`
#' - `bracket`: final bracketing interval
#' - `iters`: number of bisection iterations used
#' - `B_cal`, `alpha`: echoed inputs
#' - `trace`: data frame with the search history (gamma tried and type-I estimates)
#' - `seeds`: vector of seeds actually used (from the final evaluation)
#' - `settings`: echo of key inputs (design, prior, simulation controls)
#'
#' @examples
#' \donttest{
#' # Basic calibration without progress bar
#' cal <- bcts_calibrate_betaBinom_conj(
#'   alpha = 0.10, p_c = 0.85, M = -0.20,
#'   n_c = 29, n_t = 29,
#'   prior = "flat",
#'   B_cal = 500, n_draws = 500, seed = 11
#' )
#' cal$gamma
#' cal$type1
#' cal$trace
#'
#' # Example with a custom progress function (console)
#' pb_fun <- local({
#'   pb <- utils::txtProgressBar(min = 0, max = 10, style = 3)
#'   function(iter, maxit) {
#'     utils::setTxtProgressBar(pb, iter)
#'     if (iter == maxit) close(pb)
#'   }
#' })
#'
#' cal2 <- bcts_calibrate_betaBinom_conj(
#'   alpha = 0.10, p_c = 0.85, M = -0.20,
#'   n_c = 29, n_t = 29,
#'   prior = "flat",
#'   B_cal = 200, n_draws = 200, seed = 22,
#'   maxit = 10,
#'   show_progress = FALSE,     # disable internal bar
#'   progress_fun = pb_fun      # external progress updater
#' )
#' cal$gamma
#' cal$type1
#' cal$mc_ci
#' }
#'
#' @seealso [bcts_type1_betaBinom_conj()], [bcts_power_betaBinom_conj()]
#' @family beta-binomial simulation
#' @export
bcts_calibrate_betaBinom_conj <- function(alpha = 0.10, p_c, M, n_c, n_t,
                                          prior = c("flat","power"),
                                          prior_args = list(),
                                          B_cal = 2000, lower = 0.80, upper = 0.999,
                                          n_draws = 2000,
                                          tol = 0.001, maxit = 30,
                                          seed = 11L, show_progress = TRUE,
                                          verbose = TRUE, digits = 4,
                                          progress_fun = NULL,
                                          calibrate_on = c("point","upper","lower")) {

  prior <- match.arg(prior)
  calibrate_on <- match.arg(calibrate_on)

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0,1).")
  if (lower <= 0 || upper >= 1 || lower >= upper)
    stop("`lower` and `upper` must satisfy 0 < lower < upper < 1.")
  if (p_c < 0 || p_c > 1) stop("`p_c` must be in [0,1].")
  if (n_c < 1 || n_t < 1) stop("`n_c` and `n_t` must be positive integers.")

  set.seed(seed)
  seeds <- sample.int(1e9, B_cal)  # common random numbers

  if (is.null(tol)) {
    tol <- 1.25 * sqrt(alpha * (1 - alpha) / B_cal)
    if (verbose) message(sprintf("Auto tol set to %.4f based on B_cal=%d", tol, B_cal))
  }

  # cache full results per gamma (so we can reuse MC error, etc.)
  cache <- new.env(parent = emptyenv())
  keyfun <- function(thr) paste0("thr_", format(round(thr, digits), nsmall = digits))

  # returns the FULL list from type1 (estimate, ci_lower, ci_upper, mc_se, ...)
  type1_at <- function(thr) {
    key <- keyfun(thr)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    res <- bcts_type1_betaBinom_conj(
      B = B_cal, p_c = p_c, M = M, n_c = n_c, n_t = n_t,
      threshold = thr, prior = prior, prior_args = prior_args,
      n_draws = n_draws, seeds = seeds,
      show_progress = FALSE  # avoid nested/console bars
    )
    assign(key, res, envir = cache)
    res
  }

  #  selector for the calibration target
  get_metric <- function(res) {
    switch(calibrate_on,
           point = res$estimate,
           upper = res$ci_upper,
           lower = res$ci_lower
    )
  }

  # initial bracket diagnostics (must bracket alpha on the chosen metric)
  r_lo <- type1_at(lower); m_lo <- get_metric(r_lo)
  r_hi <- type1_at(upper); m_hi <- get_metric(r_hi)
  if (verbose) {
    message(sprintf(
      "Init: lower=%.3f -> %s=%.4f | upper=%.3f -> %s=%.4f | target alpha=%.4f",
      lower, calibrate_on, m_lo, upper, calibrate_on, m_hi, alpha))
  }
  if (!(m_lo >= alpha && m_hi <= alpha)) {
    warning(sprintf(
      "Alpha not bracketed for '%s': metric(lower)=%.4f, metric(upper)=%.4f",
      calibrate_on, m_lo, m_hi))
  }


  lo <- lower; hi <- upper
  trace <- data.frame(iter=integer(0), gamma_try=numeric(0), type1=numeric(0),
                      ci_lower = numeric(0), ci_upper = numeric(0),
                      metric=numeric(0),lo=numeric(0), hi=numeric(0), diff=numeric(0))

  # console pb only if not using Shiny progress_fun
  if (show_progress && is.null(progress_fun)) {
    pb <- utils::txtProgressBar(min = 0, max = maxit, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  if (!is.null(progress_fun)) progress_fun(0L, maxit)


  # bisection
  for (i in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    r_mid <- type1_at(mid)
    m_mid <- get_metric(r_mid)           # <- metric we calibrate on
    diff  <- m_mid - alpha

    trace <- rbind(trace, data.frame(
      iter = i, gamma_try = mid,
      type1 = r_mid$estimate,
      ci_lower  = r_mid$ci_lower,
      ci_upper  = r_mid$ci_upper,
      metric = m_mid,
      lo = lo, hi = hi, diff = diff
    ))

    if (!is.null(progress_fun)) progress_fun(i, maxit)
    if (show_progress && is.null(progress_fun)) utils::setTxtProgressBar(pb, i)
    if (verbose) message(sprintf(
      "Iter %02d: gamma=%.4f -> %s=%.4f (diff=%.4f)",
      i, mid, calibrate_on, m_mid, diff))

    if (abs(diff) < tol) {
      if (!is.null(progress_fun)) progress_fun(maxit, maxit)  # <- force 100%
      return(list(
        gamma   = mid,
        type1   = r_mid$estimate,
        mc_se   = r_mid$mc_se,
        mc_ci   = r_mid$mc_ci,
        bracket = c(lower = lo, upper = hi),
        iters   = i,
        B_cal   = B_cal,
        alpha   = alpha,
        trace   = trace,
        seeds   = r_mid$seeds,
        settings = list(
          p_c = p_c, M = M, n_c = n_c, n_t = n_t,
          prior = prior, prior_args = prior_args,
          n_draws = n_draws, tol = tol, maxit = maxit,
          seed = seed, calibrate_on = calibrate_on
        )
      ))
    }
    if (diff > 0) lo <- mid else hi <- mid
  }

  # after the loop (maxit fallback), before returning
  if (!is.null(progress_fun)) progress_fun(maxit, maxit)

  # maxit fallback
  mid  <- (lo + hi) / 2
  r_mid <- type1_at(mid)
  m_mid <- get_metric(r_mid)

  trace <- rbind(trace, data.frame(
    iter = maxit + 1, gamma_try = mid,
    type1 = r_mid$estimate, ci_lower = r_mid$ci_lower, ci_upper = r_mid$ci_upper,
    metric = m_mid, lo = lo, hi = hi, diff = m_mid - alpha))

  list(
    gamma   = mid,
    type1   = r_mid$estimate,
    mc_se   = r_mid$mc_se,
    mc_ci   = r_mid$mc_ci,
    bracket = c(lower = lo, upper = hi),
    iters   = maxit,
    B_cal   = B_cal,
    alpha   = alpha,
    trace   = trace,
    seeds   = r_mid$seeds,
    settings = list(
      p_c = p_c, M = M, n_c = n_c, n_t = n_t,
      prior = prior, prior_args = prior_args,
      n_draws = n_draws, tol = tol, maxit = maxit,
      seed = seed, calibrate_on = calibrate_on
    )
  )
}
