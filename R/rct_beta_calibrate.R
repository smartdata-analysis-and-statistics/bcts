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
#' @param prior Character. Specifies the prior distribution. Options:
#'   - `"flat"`: Flat Beta(1,1) prior for both arms (non-informative),
#'   - `"jeffreys"`: Jeffreys prior, i.e., Beta(0.5, 0.5) for both arms,
#'   - `"power"`: Power prior on the control arm, requires specification in `prior_args`.
#' @param prior_args List of additional prior hyperparameters. Required when `prior = "power"`:
#'   - `a0`: Discount factor for historical control data (numeric in `[0,1]`),
#'   - `y_0`, `n_0`: Historical control responders and total sample size (integers),
#'   - `a_base`, `b_base`: Baseline Beta parameters (default is `1`, `1`).
#'   Ignored when `prior` is `"flat"` or `"jeffreys"`.
#' @param B_cal Number of simulated trials for calibration of `gamma`.
#' @param B_power Number of simulated trials for power estimation.
#' @param n_draws Posterior Monte Carlo draws per trial for evaluating `Pr(NI)`.
#' @param tol Absolute tolerance for calibration convergence.
#' @param seed RNG seed for reproducibility.
#' @param show_progress Logical. Show text progress bars during simulation.
#' @param verbose Logical. Print iteration messages during calibration.
#' @param method Character. Method for estimation. Either `"simulate"` or `"cpp"`.
#'
#' @return An object of class `"bayesNI"` with components:
#' - `calibration`: list from [bcts_calibrate_betaBinom_conj()]
#' - `power`: estimated Bayesian power at calibrated `gamma`
#' - `settings`: list of all input parameters (echoed)
#'
#' @seealso [bcts_calibrate_betaBinom_conj()],
#'   [rct_beta_power()],
#'   [bcts_type1_betaBinom_conj()]
#' @family conjugate-NI
#'
#' @examples
#' \donttest{
#' res <- rct_beta_calibrate(
#'   alpha = 0.10, p_c = .85, p_t = .85,
#'   n_c = 29, n_t = 29, NI_margin = -0.20,
#'   prior = "flat", B_cal = 1000, B_power = 500,
#'   method = "cpp"
#' )
#' res$calibration$gamma
#' res$power
#' }
#'
#' @export
rct_beta_calibrate <- function(alpha = 0.10,
                              p_c = 0.85, p_t = 0.85,
                              n_c = 29, n_t = 29,
                              NI_margin = -0.20,
                              prior = c("flat", "jeffreys", "power"),
                              prior_args = list(),
                              B_cal = 2000, B_power = 1000,
                              n_draws = 2000,
                              tol = 0.002,
                              seed = 123,
                              show_progress = FALSE, verbose = TRUE,
                              method = c("cpp", "simulate")) {
  prior <- match.arg(prior)

  method <- match.arg(method)

  # 1) Calibrate gamma using conjugate machinery
  cal <- bcts_calibrate_betaBinom_conj(alpha = alpha, p_c = p_c, M = NI_margin,
                              n_c = n_c, n_t = n_t,
                              prior = prior, prior_args = prior_args,
                              B_cal = B_cal, n_draws = n_draws,
                              tol = tol, seed = seed,
                              show_progress = show_progress, verbose = verbose)

  # 2) Simulate Bayesian "power" at calibrated gamma
  power <- rct_beta_power(B = B_power, p_c = p_c, p_t = p_t,
                          n_c = n_c, n_t = n_t, M = NI_margin,
                          threshold = cal$gamma,
                          prior = prior, prior_args = prior_args,
                          n_draws = n_draws,
                          show_progress = show_progress,
                          method = method, seed = seed)


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


#' Plot calibration trace for calibrated Bayesian NI threshold
#'
#' Visualizes the estimated Type-I error at different posterior thresholds
#' during the calibration procedure.
#'
#' @param x Output from [bcts_calibrate_betaBinom_conj()]
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#'
#' @examples
#' # res <- bcts_calibrate_betaBinom_conj(...)  # Run separately
#' # plot(res)
plot.bcts_calibration <- function(x, ...) {
  stopifnot(!is.null(x$trace), !is.null(x$type1))

  tr     <- x$trace
  alpha  <- x$alpha
  gamma  <- x$gamma
  type1  <- x$type1

  # Small vertical offset above CI
  offset <- 0.01

  ggplot2::ggplot(tr, ggplot2::aes(x = .data$gamma_try, y = .data$type1)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      width = 0.002, alpha = 0.4
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      x = gamma,
      y = type1["estimate"],
      color = "blue",
      size = 3
    ) +
    ggplot2::annotate(
      "text",
      x = gamma,
      y = type1["ci_upper"] + offset,
      label = sprintf("gamma = %.3f\nType-I = %.1f%%",
                      gamma, 100 * type1["estimate"]),
      hjust = 0.5,
      size = 3.5
    ) +
    ggplot2::geom_hline(
      yintercept = alpha,
      linetype = "dotted",
      color = "red"
    ) +
    ggplot2::labs(
      x = expression(gamma),
      y = "Estimated Type-I error",
      subtitle = sprintf(
        "Calibrated gamma = %.3f after %d steps; dotted = target alpha = %.2f",
        gamma, x$iters, alpha
      )
    ) +
    ggplot2::theme_minimal(base_size = 12)+
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1))+
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1))
}
