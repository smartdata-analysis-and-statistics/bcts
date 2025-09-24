#' Estimate Bayesian Power for a 2-arm Binomial RCT (Beta-Binomial Conjugate)
#'
#' @description Computes the power of a Bayesian randomized controlled trial
#' using conjugate Beta-Binomial models. Can use exact simulation in R or fast
#' simulation in C++.
#'
#' @param B Integer. Number of trial simulations (ignored if \code{method = "exact"}).
#' @param p_c,p_t Numeric in \[0,1\]. True response probabilities in control and treatment arms.
#' @param n_c,n_t Integers. Sample sizes in control and treatment arms.
#' @param M Numeric. Margin on the risk-difference scale. Default is 0 (superiority/equivalence).
#' @param threshold Posterior probability cutoff (e.g. 0.9).
#' @param prior Character. Prior type: \code{"flat"} or \code{"power"}.
#' @param prior_args List with prior hyperparameters. Optional elements:
#'   - \code{a0}, \code{y_0}, \code{n_0} for power prior
#'   - \code{a_base_c}, \code{b_base_c} for control arm prior
#'   - \code{a_base_t}, \code{b_base_t} for treatment arm prior
#' @param n_draws Integer. Number of posterior samples per trial.
#' @param show_progress Logical. Show text progress bar?
#' @param method Character. Either \code{"simulate"} (pure R) or \code{"cpp"} (fast C++).
#' @param seed Optional integer. If provided, sets the RNG seed for reproducibility.
#'
#' @return A list with:
#' \item{estimate}{Estimated power (probability of success)}
#' \item{mc_se}{Monte Carlo standard error}
#' \item{successes}{Number of trials with successful outcome}
#' \item{B}{Number of simulations}
#'
#' @examples
#' # Using C++ backend
#' rct_beta_power(
#'   B = 1000, p_c = 0.85, p_t = 0.85,
#'   n_c = 30, n_t = 30, M = -0.1,
#'   threshold = 0.9, prior = "flat",
#'   method = "cpp", show_progress = FALSE
#' )
#'
#' @export
rct_beta_power <- function(B = 10000,
                           p_c, p_t, n_c, n_t,
                           M = 0,
                           threshold,
                           prior = c("flat", "power"),
                           prior_args = list(),
                           n_draws = 2000,
                           show_progress = TRUE,
                           method = c("cpp", "simulate"),
                           seed = NULL) {
  prior <- match.arg(prior)
  method <- match.arg(method)

  if (!is.null(seed)) set.seed(seed)

  if (method == "simulate") {
    res <- bcts_power_betaBinom_conj(
      B = B,
      p_c = p_c,
      p_t = p_t,
      n_c = n_c,
      n_t = n_t,
      M = M,
      threshold = threshold,
      prior = prior,
      prior_args = prior_args,
      n_draws = n_draws,
      show_progress = show_progress,
      seed = seed
    )
    return(res)
  }

  if (method == "cpp") {
    # Ensure required C++ function is available
    dec_vec <- rct_power_beta_binom_cpp_vec(
      B = B,
      p_c = p_c,
      p_t = p_t,
      n_c = n_c,
      n_t = n_t,
      M = M,
      threshold = threshold,
      prior = prior,
      prior_args = prior_args,
      n_draws = n_draws,
      show_progress = show_progress
    )
    k <- sum(dec_vec)
    p <- mean(dec_vec)
    mc_se <- sqrt(p * (1 - p) / B)

    ci <- binom::binom.confint(x=k, n=B, conf.level = 0.95, methods = "wilson")

    return(list(
      estimate = p,
      ci_lower = ci["lower"],
      ci_upper = ci["upper"],
      mc_se = mc_se,
      successes = k,
      B = B
    ))
  }

  stop("Unsupported method. Choose either 'simulate' or 'cpp'.")
}
