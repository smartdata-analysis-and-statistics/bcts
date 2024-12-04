#' Simulate a randomized clinical trial
#'
#' @param n Sample size
#' @param mean vector with the means of the treatment groups
#' @param sd vector with the standard deviations of the treatment groups
#' @param trtnames vector with treatment names
#'
#' @return A dataframe with the simulated trial data
#'
#' @author Thomas Debray
#'
#' @export
#' @importFrom stats rnorm
#' @importFrom dplyr mutate rename slice_head
#' @importFrom rlang .data
#' @importFrom randomizeR crPar genSeq getRandList
#'
sim_rct_normal <- function(n,
                           mean,
                           sd,
                           trtnames) {

  if (missing(trtnames)) {
    if (is.null(names(mean))) {
      stop("No treatment names provided")
    }
    trtnames = names(mean)
  }

  # Generate randomization parameters
  params <- crPar(N = n,
                  K = length(trtnames),
                  ratio = rep(1, length(trtnames)),
                  groups = trtnames)

  # Get randomization sequence
  rs <- genSeq(params)

  dat <- data.frame(id = seq(n),
                    Treatment = getRandList(rs)[1,],
                    Y = NA
                    )

  if (any(is.na(c(mean,sd)))) {
    dat <- dat %>% mutate(Treatment = factor(.data$Treatment))
    return(dat)
  }


  # Create a data frame to store the simulated trial data
  for (trt in seq(trtnames)) {
    ntrt <- sum(dat$Treatment == trtnames[trt])
    dat$Y[which(dat$Treatment == trtnames[trt])] <- rnorm(n = ntrt, mean = mean[trt], sd[trt])
  }
  dat <- dat %>% mutate(Treatment = factor(.data$Treatment))

  return(dat)
}

#' Monte Carlo Error for Proportion Estimate
#'
#' This function calculates the Monte Carlo error for a proportion estimate based on observed event counts.
#'
#' @param x Numeric. The observed event count.
#' @param n Numeric. The total number of Monte Carlo simulations.
#' @param level Numeric. The confidence level for the error bounds (default is 0.95).
#' @return A data frame containing the observed event count, total simulations, proportion estimate, standard error, and confidence interval.
#' @importFrom rlang .data
#' @export
#'
#' @author Thomas Debray
mc_error_proportion <- function(x, n, level = 0.95) {
  out <- data.frame(x = x,
                    n = n,
                    est = x/n) %>% dplyr::mutate(
    se = sqrt((.data$est*(1 - .data$est))/n),
    lower = .data$est + stats::qnorm((1 - level)/2)*.data$se,
    upper = .data$est + stats::qnorm( 1 - (1 - level)/2)*.data$se)

  return(out)
}

#' @title Gamma Calibration for Target Type-I Error
#' @description
#' This function calibrates the gamma threshold in a Bayesian clinical trial simulation to achieve a specified target Type-I error rate.
#' It uses a binomial regression model to iteratively refine gamma, aiming for the closest match to the target Type-I error within a given tolerance.
#'
#' @param bcts_fun A function that simulates Bayesian clinical trials and computes results for a given gamma threshold.
#'                 It should return a list containing a `simresults` object with details such as `gamma`, `nsim`, and `rejectH0.final`.
#' @param type1 The target Type-I error rate.
#' @param type1.tolerance The acceptable tolerance for the difference between the empirical and target Type-I error rate. Default is `0.001`.
#' @param nsim Number of trials to simulate (default: 1000)
#' @param gamma_scaling A numeric value used to scale or adjust the gamma parameter during the initial step of the Type-I error estimation process.
#'                      This factor is applied as an exponent to `1-gamma` in initial calculations.
#' @param req_n_events Initial number of required events for the simulation. This value is increased iteratively to ensure robust calibration. Default is `10`.
#' @param max.iter The maximum number of iterations allowed during the calibration process. Default is `10`.
#'
#' @details
#' The function begins by estimating the Type-I error for `(1-type1)` and `(1-type1)^gamma_scaling`.
#' It then fits a binomial regression model to predict the gamma threshold that achieves the target Type-I error.
#' During each iteration, the function:
#' - Computes the empirical Type-I error for the predicted gamma threshold.
#' - Updates the binomial regression model with the new results.
#' - Refines the gamma threshold until the empirical Type-I error is within the specified tolerance or the maximum number of iterations is reached.
#'
#' If the tolerance is not achieved within `max.iter` iterations, the function returns the best available estimate for the gamma threshold.
#'
#' @return
#' A list containing the following elements:
#' - **gamma.opt**: The calibrated gamma threshold for declaring success, optimized to approximate the target Type-I error.
#' - **sim.opt**: The simulation results for the final calibrated gamma threshold.
#' - **table**: A `data.frame` summarizing all tested gamma threshold values and their corresponding results:
#'     - `gamma`: The tested gamma threshold values.
#'     - `n.exp`: The expected number of events based on the tested gamma threshold.
#'     - `n.obs`: The observed number of events.
#'     - `n.total`: The total number of simulations.
#'     - `est`: The estimated Type-I error rate.
#'     - `lower`: The lower bound of the confidence interval for the Type-I error.
#'     - `upper`: The upper bound of the confidence interval for the Type-I error.
#'
#' @note
#' This function assumes that `bcts_fun` returns a list containing a `simresults` object with the following:
#' - `gamma`: The gamma threshold value used in the simulation.
#' - `nsim`: The number of simulations performed.
#' - `rejectH0.final`: A vector indicating whether the null hypothesis was rejected in each simulation.
#'
#' The binomial regression model is used to interpolate the optimal gamma threshold rather than relying solely on stepwise iteration, improving efficiency and precision.
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @importFrom stats glm binomial coef
#' @importFrom dplyr add_row %>%

calibrate_gamma <- function(bcts_fun,
                            type1 = 0.025,
                            type1.tolerance = 0.001,
                            nsim = 1000,
                            gamma_scaling = 0.75,
                            req_n_events = 10,
                            max.iter = 10) {
  iter <- 1

  # Estimate alpha for the initial gamma
  sim.new <- bcts_fun(1 - type1, nsim = ceiling(req_n_events/type1))
  power.sim.new <- power(sim.new, adjust_for_futility = FALSE) #No futility adjustment for assessing type-1

  # Estimate alpha for an extreme value of gamma
  sim.cal <- bcts_fun((1 - type1)^gamma_scaling, nsim = ceiling(req_n_events/(1 - (1 - type1)^gamma_scaling)))
  power.sim.cal <- power(sim.cal, adjust_for_futility = FALSE) #No futility adjustment for assessing type-1

  # Estimate binomial model
  ds.bin <- data.frame(gamma = c(sim.cal$gamma, sim.new$gamma),
                       n.exp = c((1 - sim.cal$gamma)*sim.cal$nsim, (1 - sim.new$gamma)*sim.new$nsim),
                       n.obs = c(power.sim.cal$x, power.sim.new$x),
                       n.total = c(power.sim.cal$n, power.sim.new$n),
                       est = c(power.sim.cal$est, power.sim.new$est),
                       lower = c(power.sim.cal$lower, power.sim.new$lower),
                       upper = c(power.sim.cal$upper, power.sim.new$upper))

  while (abs(power.sim.new$est - type1) > type1.tolerance & iter < max.iter) {

    fit <- glm(cbind(n.obs, n.total-n.obs) ~ log(gamma/(1-gamma)),
               family = binomial(), data = ds.bin)

    # exp(beta0 + beta1* Tgamma) = p.success
    # beta0 + beta1* Tgamma = log(p.success)
    # Tgamma = (log(p.success)-beta0)/beta1

    gamma.test.new <- 1/(1 + exp(-(log(type1) - coef(fit)[1])/coef(fit)[2]))
    sim.new <- bcts_fun(gamma.test.new, nsim = nsim)
    power.sim.new <- power(sim.new, adjust_for_futility = FALSE)

    ds.bin <- ds.bin %>% add_row(
      data.frame(gamma = sim.new$gamma,
                 n.exp = (1 - sim.new$gamma)*sim.new$nsim,
                 n.obs = power.sim.new$x,
                 n.total = power.sim.new$n,
                 est = power.sim.new$est,
                 lower = power.sim.new$lower,
                 upper = power.sim.new$upper)
      )

    iter <- iter + 1

    # Increase required number of events
    req_n_events <- req_n_events*2.5
  }

  if (iter == max.iter) {
    warning("Maximum number of iterations reached. The calibration may not have converged.")
  }

  return(list(gamma.opt = sim.new$gamma,
              sim.opt = sim.new,
              table = ds.bin))
}

