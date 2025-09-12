# Draw-based approximation to Pr(theta_t - theta_c > M)
# given posterior Beta shapes.
rbeta_diff_prob <- function(a_t, b_t, a_c, b_c, M, n_draws = 2000) {
  pt <- stats::rbeta(n_draws, a_t, b_t)
  pc <- stats::rbeta(n_draws, a_c, b_c)
  mean((pt - pc) > M)
}
# Build posterior Beta shapes for a single dataset under chosen prior.
# Power prior applies to CONTROL only (as in your JAGS version).
posterior_shapes_conj <- function(y_c, n_c, y_t, n_t,
                                  prior = c("flat","power"),
                                  prior_args = list()) {
  prior <- match.arg(prior)

  if (prior == "flat") {
    # Beta(1,1) prior for both arms
    a_c <- 1 + y_c
    b_c <- 1 + (n_c - y_c)
    a_t <- 1 + y_t
    b_t <- 1 + (n_t - y_t)
  } else {
    a0     <- prior_args$a0     %||% 0.5
    y_0    <- prior_args$y_0    %||% stop("prior_args$y_0 required for power prior")
    n_0    <- prior_args$n_0    %||% stop("prior_args$n_0 required for power prior")
    a_base <- prior_args$a_base %||% 1
    b_base <- prior_args$b_base %||% 1

    # Control arm: power prior Beta(a_base + a0*y0, b_base + a0*(n0-y0))
    a_c <- a_base + a0 * y_0 + y_c
    b_c <- b_base + a0 * (n_0 - y_0) + (n_c - y_c)

    # Treatment arm: flat Beta(1,1)
    a_t <- 1 + y_t
    b_t <- 1 + (n_t - y_t)
  }

  list(a_c = a_c, b_c = b_c, a_t = a_t, b_t = b_t)
}

#' Simulate one NI trial (Beta–Binomial, conjugate)
#'
#' Simulates one two-arm non-inferiority trial with binomial outcomes, updates the
#' Beta posteriors by conjugacy (no MCMC), and returns the posterior probability
#' of non-inferiority \eqn{Pr(\theta_t - \theta_c > M \mid y)}.
#'
#' @param p_c,p_t Numeric in \[0,1\]. True response probabilities for control and treatment.
#' @param n_c,n_t Integers. Sample sizes for control and treatment.
#' @param M Numeric. NI margin on the risk-difference scale (e.g., `-0.20`).
#' @param prior Character. `"flat"` (Beta(1,1) both arms) or `"power"` (power prior on control).
#' @param prior_args List of hyperparameters when `prior = "power"`:
#'   - `a0`: discount factor in \[0,1\]
#'   - `y_0`, `n_0`: historical control responders and sample size
#'   - `a_base`, `b_base`: baseline Beta parameters (defaults `1`, `1`)
#' @param n_draws Integer. Monte Carlo draws used to evaluate `Pr(NI)` from the two Beta posteriors.
#' @param seed Optional integer RNG seed for this single trial.
#'
#' @return A list with components:
#' - `summary`: named numeric vector with
#'   - `post_prob_NI`: posterior probability of NI
#'   - `mean_theta_c`, `mean_theta_t`: posterior means of control/treatment
#'   - `y_c`, `y_t`: simulated successes
#' - `shapes`: list of Beta shapes `a_c,b_c,a_t,b_t` used for the posteriors
#'
#' @seealso [bcts_power_betaBinom_conj()], [bcts_type1_betaBinom_conj()],
#'   [bcts_calibrate_betaBinom_conj()]
#'
#' @examples
#' # \donttest{
#' set.seed(123)
#' out <- bayesNI_trial_betaBinom_conj(
#'   p_c = .85, p_t = .85, n_c = 29, n_t = 29, M = -0.20,
#'   prior = "flat", n_draws = 2000
#' )
#' out$summary["post_prob_NI"]
#' # }
#'
#' @export
bayesNI_trial_betaBinom_conj <- function(p_c, p_t, n_c, n_t, M,
                                         prior = c("flat","power"),
                                         prior_args = list(),
                                         n_draws = 2000,
                                         seed = NULL) {
  prior <- match.arg(prior)
  if (!is.null(seed)) set.seed(seed)

  # Simulate trial data
  y_c <- stats::rbinom(1L, n_c, p_c)
  y_t <- stats::rbinom(1L, n_t, p_t)

  # Posterior Beta shapes for both arms under chosen prior
  sh <- posterior_shapes_conj(
    y_c = y_c, n_c = n_c, y_t = y_t, n_t = n_t,
    prior = prior, prior_args = prior_args
  )

  # Compute Pr(theta_t - theta_c > M) by MC from the two Betas
  post_prob_NI <- rbeta_diff_prob(
    a_t = sh$a_t, b_t = sh$b_t,
    a_c = sh$a_c, b_c = sh$b_c,
    M = M, n_draws = n_draws
  )

  est <- c(
    post_prob_NI = post_prob_NI,
    mean_theta_c = sh$a_c / (sh$a_c + sh$b_c),
    mean_theta_t = sh$a_t / (sh$a_t + sh$b_t),
    y_c = y_c, y_t = y_t
  )

  list(summary = est, shapes = sh)
}
