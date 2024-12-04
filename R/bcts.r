#' Bayesian Clinical Trial Simulation (BCTS)
#'
#' @param n_dose_sel First interim sample size (for dose selection)
#' @param n_ss_reest Optional second interim sample size (for sample size re-estimation)
#' @param n_pln Planned sample size
#' @param n_max Maximum sample size
#' @param mu Named vector with expected effect sizes
#' @param sigma Named vector with expected standard deviations
#' @param trt_ref Character denoting the control treatment
#' @param trt_rank Named vector with preference ranking for each treatment (e.g., lower doses are preferred) when multiple treatments have a posterior predictive power > 'th.eff'
#' @param alpha Desired one-sided Type-I error rate. Default is `0.025`.
#'              This represents the maximum allowable probability of rejecting the null hypothesis when it is true.
#' @param alpha_tolerance Tolerance for the Type-I error rate. Default is `0.001`.
#'                        This specifies the acceptable deviation from the target `alpha` during the calibration process.
#'                        Smaller values indicate stricter adherence to the target, while larger values allow more flexibility.
#' @param th.fut Futility threshold for predictive power (default: 0.2)
#' @param th.eff Efficacy threshold for predictive power (default: 0.9)
#' @param th.prom Predictive power threshold to trigger sample size increase (default: 0.5)
#' @param method Method to estimate the power ("bayes" or "mcmc", default: "mcmc")
#' @param nsim Number of trials to simulate (default: 1000)
#' @param num_chains Number of MCMC chains (default: 4)
#' @param n.iter Total MCMC iterations (default: 5000)
#' @param n.adapt MCMC adaptation iterations (default: 500)
#' @param perc_burnin Proportion of iterations for burn-in (default: 0.2)
#' @param progress.bar Progress bar type ("text", "gui", "none")
#'
#' @description
#' Simulates Bayesian clinical trials with up to two interim analyses for dose selection
#' and sample size re-estimation. The null hypothesis is rejected if the posterior probability
#' $Pr(\\mu_t - \\mu_c > 0 | D) > \\gamma$ exceeds a defined threshold.
#'
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @return An object of class "bcts"
#' @export
bcts <- function(n_dose_sel, n_ss_reest, n_pln, n_max, mu, sigma,
                 trt_ref,
                 trt_rank = NULL, alpha = 0.025, alpha_tolerance = 0.001,
                 th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                 method = "mcmc", nsim = 1000, num_chains = 4,
                 n.iter = 5000, n.adapt = 500, perc_burnin = 0.2,
                 progress.bar = "text") {


  # Input validation
  stopifnot(is.numeric(n_dose_sel), is.numeric(n_ss_reest),
            is.numeric(n_pln), is.numeric(n_max),
            is.numeric(alpha), is.numeric(th.fut), is.numeric(th.eff),
            is.numeric(th.prom), method %in% c("bayes", "mcmc"))

  # Set effect sizes for evaluating type-I error
  mu.type1 <- mu * 0

  if (is.null(n_ss_reest)) {
    # Calibrate gamma for a one-interim design
    message("Using bcts_one_interim.")
    opt.gamma <- calibrate_gamma(
      bcts_fun = function(x, req_n_events) {
        bcts_one_interim(n_int = n_dose_sel, n_pln = n_pln, n_max = n_max,
                         mu = mu.type1, sigma = sigma, trt_ref = trt_ref,
                         gamma = x,
                         trt_rank = trt_rank, th.fut = th.fut, th.eff = th.eff,
                         th.prom = th.prom, method = method,
                         nsim = ceiling(req_n_events/(1 - x)),
                         num_chains = num_chains, n.iter = n.iter,
                         n.adapt = n.adapt, perc_burnin = perc_burnin,
                         progress.bar = progress.bar)
      }, type1 = alpha, type1.tolerance = alpha_tolerance)

    ## Assess type-1 error
    sim_type1 <- bcts_one_interim(n_int = n_dose_sel, n_pln = n_pln,
                                  n_max = n_max,
                                  mu = mu.type1, sigma = sigma,
                                  trt_ref = trt_ref,
                                  gamma = opt.gamma$gamma.opt,
                                  trt_rank = trt_rank, th.fut = th.fut, th.eff = th.eff,
                                  th.prom = th.prom, method = method, nsim = nsim,
                                  num_chains = num_chains, n.iter = n.iter,
                                  n.adapt = n.adapt, perc_burnin = perc_burnin,
                                  progress.bar = progress.bar)

    ## Assess power
    sim_power <- bcts_one_interim(n_int = n_dose_sel, n_pln = n_pln,
                                  n_max = n_max,
                                  mu = mu, sigma = sigma, trt_ref = trt_ref,
                                  gamma = opt.gamma$gamma.opt,
                                  trt_rank = trt_rank, th.fut = th.fut, th.eff = th.eff,
                                  th.prom = th.prom, method = method, nsim = nsim,
                                  num_chains = num_chains, n.iter = n.iter,
                                  n.adapt = n.adapt, perc_burnin = perc_burnin,
                                  progress.bar = progress.bar)
  } else if (!is.null(n_ss_reest) & !is.null(n_dose_sel)) {
    # Calibrate gamma for the two-interim version
    message("Using bcts_two_interim.")
    opt.gamma <- calibrate_gamma(
      bcts_fun = function(x, req_n_events) {
        bcts_two_interim(n_int1 = n_dose_sel, n_int2 = n_ss_reest,
                         n_pln = n_pln, n_max = n_max,
                         mu = mu.type1, sigma = sigma, trt_ref = trt_ref,
                         gamma = x, trt_rank = trt_rank, th.fut = th.fut,
                         th.eff = th.eff, th.prom = th.prom, method = method,
                         nsim = ceiling(req_n_events/(1 - x)),
                         num_chains = num_chains, n.iter = n.iter,
                         n.adapt = n.adapt, perc_burnin = perc_burnin,
                         progress.bar = progress.bar)
      }, type1 = alpha, type1.tolerance = alpha_tolerance)

    ## Assess type-1 error
    sim_type1 <- bcts_two_interim(n_int1 = n_dose_sel, n_int2 = n_ss_reest,
                                  n_pln = n_pln, n_max = n_max,
                                  mu = mu.type1, sigma = sigma, trt_ref = trt_ref,
                                  gamma = opt.gamma$gamma.opt,
                                  trt_rank = trt_rank, th.fut = th.fut,
                                  th.eff = th.eff, th.prom = th.prom, method = method,
                                  nsim = nsim,
                                  num_chains = num_chains, n.iter = n.iter,
                                  n.adapt = n.adapt, perc_burnin = perc_burnin,
                                  progress.bar = progress.bar)

    ## Assess power
    sim_power <- bcts_two_interim(n_int1 = n_dose_sel, n_int2 = n_ss_reest,
                                  n_pln = n_pln, n_max = n_max,
                                  mu = mu, sigma = sigma, trt_ref = trt_ref,
                                  gamma = opt.gamma$gamma.opt,
                                  trt_rank = trt_rank, th.fut = th.fut,
                                  th.eff = th.eff, th.prom = th.prom, method = method,
                                  nsim = nsim,
                                  num_chains = num_chains, n.iter = n.iter,
                                  n.adapt = n.adapt, perc_burnin = perc_burnin,
                                  progress.bar = progress.bar)

  } else {
    stop("Please provide either `n_ss_reest` or `n_dose_sel`.")
  }


  out <- list(sim_type1 = sim_type1, sim_power = sim_power,
              alpha = alpha,
              opt.gamma = opt.gamma)
  class(out) <- "bcts_results"

  return(out)
}



#' Estimate the Power for the Difference between Two Means
#' @param n1 Number of observations in treatment group 1
#' @param n2 Number of observations in treatment group 2
#' @param N Total sample size at final analysis (for calculating PPoS)
#' @param mu1 Expected mean in treatment group 1
#' @param mu2 Expected mean in treatment group 2
#' @param margin Margin for deciding superiority/non-inferiority
#' @param sd1 Standard deviation of the outcome in treatment group 1
#' @param sd2 Standard deviation of the outcome in treatment group 2
#' @param alpha One-sided alpha
#' @param alternative Alternative hypothesis, defaults to "Superiority"
#'
#' @description
#' For superiority, H0 is defined as mu1 - mu2 <= margin.
#'
#' @importFrom stats pnorm qnorm
#'
#' @return An object of class "fpower"
#'
pow_diff_means <- function(n1, n2, N, mu1, mu2, margin = 0,
                             sd1,
                             sd2,
                             alpha,
                             alternative = "superiority") {

  # Calculate difference
  mudiff <- mu1 - mu2 - margin

  # calculate pooled SD
  sd_pooled <- sqrt(((n1 - 1)*sd1**2 + (n2 - 1)*sd2**2)/(n1 + n2 - 2))

  # Calculate pooled SE
  se_pooled <- sd_pooled * sqrt(1/n1 + 1/n2)

  # Calculate allocation ratio
  alr <- n1/n2
  r <- sqrt(((alr + 1)^2) / alr)

  # Calculate total sample size for which observations are available
  n <- n1 + n2

  if (alternative == "superiority") {
    z_test <- (mu1 - mu2 - margin)/se_pooled
    z_crit <- qnorm(1 - alpha)

    # Calculate the probability that mu1 - mu2 - margin > 0
    power <- pnorm(z_test - z_crit)

    # Calculate Predictive Power of success (PPoS)
    ppos <- pnorm((1/(r*sd_pooled))*sqrt(n/(N-n))* ((mu1 - mu2 - margin)*sqrt(N)-r*sd_pooled*z_crit))


    str_out <- "Difference between Two means\n"
    str_out <- paste0(str_out, "(Test for Noninferority/Superiority)\n")
    str_out <- paste0(str_out, "H0: mu1 - mu2 <= ", margin, "\n")
    str_out <- paste0(str_out, "H0: mu1 - mu2 > ", margin, "\n")
    str_out <- paste0(str_out, "------------------------------ \n")
    str_out <- paste0(str_out, " Statistical power = ", round(power,3), "\n")
    str_out <- paste0(str_out, " n1 = ", n1, "\n")
    str_out <- paste0(str_out, " n2 = ", n2, "\n")

    out <- list(diff = mudiff,
                diff_lower = mudiff + qnorm(alpha) * se_pooled,
                diff_upper = mudiff + qnorm(1 - alpha) * se_pooled,
                diff_se = se_pooled,
                power = power, ppos = ppos, Z_test = z_test, Z_crit = z_crit,
                I = 1/(se_pooled**2),
                rejectH0 = z_test > z_crit, txt = str_out)
    class(out) <- "fpower"
    return(out)
  }
}

#' Plot posterior distribution of the mean outcome
#'
#' Function to plot the posterior distribution of the mean outcome in the selected treatment arm
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export posterior
posterior <- function(...) {
  UseMethod("posterior")
}

#' Extract the power of a bcts simulation
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export power
power <- function(...) {
  UseMethod("power")
}

#' Print bcts output
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @importFrom binom binom.confint
#' @importFrom dplyr %>%
#' @importFrom rlang .data
print.bcts <- function(x, ...) {

  conf.level <- 0.95

  ## Get quantiles of final sample size
  n_fin_qt <- stats::quantile(x$simresults$n.final, c((1 - conf.level)/2, 1 - (1 - conf.level)/2))


  power_without_fut <- power(x = x, adjust_for_futility = FALSE, level = conf.level)
  power_with_fut <- power(x = x, adjust_for_futility = TRUE, level = conf.level)

  fut.trig <- mc_error_proportion(x = sum(x$simresults$fut.trig),
                                  n = nrow(x$simresults),
                                  level = conf.level)

  inc.ss <- mc_error_proportion(x = sum(x$simresults$inc.ss),
                                n = nrow(x$simresults),
                                level = conf.level)

  out <- data.frame(statistic = character(),
                      est = numeric(),
                      cil = numeric(),
                      ciu = numeric())
  out <- out %>% add_row(data.frame(statistic = "Power (Without Futility Adjustment)",
                                      est = power_without_fut$est,
                                      cil = power_without_fut$lower,
                                      ciu = power_without_fut$upper))
  out <- out %>% add_row(data.frame(statistic = "Power (With Futility Adjustment)",
                                    est = power_with_fut$est,
                                    cil = power_with_fut$lower,
                                    ciu = power_with_fut$upper))
    out <- out %>% add_row(data.frame(statistic = "Pr(fut.trig)",
                                      est = fut.trig$est,
                                      cil = fut.trig$lower,
                                      ciu = fut.trig$upper))
    out <- out %>% add_row(data.frame(statistic = "Pr(inc.ss)",
                                      est = inc.ss$est,
                                      cil = inc.ss$lower,
                                      ciu = inc.ss$upper))
    out <- out %>% add_row(data.frame(statistic = "N final",
                                      est = mean(x$simresults$n.final),
                                      cil = n_fin_qt[1],
                                      ciu = n_fin_qt[2]))

    cat(paste("Sample size calculation for an adaptive trial with", x$no.looks, "looks.\n\n"))
    print(out)
}

#' @title Print Summary of bcts_results
#' @description
#' This function provides a detailed summary of an object of class `bcts_results`, which contains simulation results for an adaptive clinical trial.
#' It computes and displays key statistics such as Type-I error, power (with and without futility adjustment), probabilities of triggering futility or increasing sample size, and final sample size estimates.
#'
#' @param x An object of class `bcts_results` containing the results of Type-I error and power simulations.
#'          The object should include components such as `sim_type1` and `sim_power`, each containing simulation results.
#' @param ... Optional additional arguments (currently unused).
#'
#' @details
#' The function extracts and summarizes important metrics from the `bcts_results` object:
#' - **Type-I Error**: The estimated Type-I error rate along with its confidence interval.
#' - **Power**:
#'   - Without futility adjustment: The power estimated without considering futility triggers.
#'   - With futility adjustment: The power adjusted for scenarios where futility was triggered.
#' - **Pr(fut.trig)**: The probability of triggering futility during the trial.
#' - **Pr(inc.ss)**: The probability of increasing the sample size during the trial.
#' - **N final**: The mean final sample size and its confidence interval.
#'
#' @return
#' A printed summary table displaying:
#' - `statistic`: The name of the statistic being summarized.
#' - `est`: The estimated value of the statistic.
#' - `cil`: The lower bound of the 95% confidence interval (or specified level).
#' - `ciu`: The upper bound of the 95% confidence interval (or specified level).
#'
#' @importFrom binom binom.confint
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom knitr kable
#' @export
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
print.bcts_results <- function(x, ...) {

  conf.level <- 0.95

  ## Get quantiles of final sample size
  n_fin_qt <- stats::quantile(x$sim_power$simresults$n.final, c((1 - conf.level)/2, 1 - (1 - conf.level)/2))

  type1 <- power(x = x$sim_type1, adjust_for_futility = FALSE, level = conf.level)

  power_without_fut <- power(x = x$sim_power, adjust_for_futility = FALSE, level = conf.level)
  power_with_fut <- power(x = x$sim_power, adjust_for_futility = TRUE, level = conf.level)

  fut.trig <- mc_error_proportion(x = sum(x$sim_power$simresults$fut.trig),
                                  n = nrow(x$sim_power$simresults),
                                  level = conf.level)

  inc.ss <- mc_error_proportion(x = sum(x$sim_power$simresults$inc.ss),
                                n = nrow(x$sim_power$simresults),
                                level = conf.level)

  cat("Sample Size Calculation for an Adaptive Trial\n")
  cat("------------------------------------------------------------\n")
  cat(paste("Target Type-I Error:", x$alpha*100, "%\n"))
  cat(paste("Number of Looks:", x$sim_type1$no.looks))

  out <- data.frame(
    Statistic = c("Type-I Error",
                  "Power (Without Futility Adjustment)",
                  "Power (With Futility Adjustment)",
                  "Pr(Futility Triggered)",
                  "Pr(Sample Size Increased)",
                  "Final Sample Size (Mean)"),
    Estimate = c(type1$est, power_without_fut$est, power_with_fut$est,
                 fut.trig$est, inc.ss$est, mean(x$sim_power$simresults$n.final)),
    `CI Lower` = c(type1$lower, power_without_fut$lower, power_with_fut$lower,
                   fut.trig$lower, inc.ss$lower, n_fin_qt[1]),
    `CI Upper` = c(type1$upper, power_without_fut$upper, power_with_fut$upper,
                   fut.trig$upper, inc.ss$upper, n_fin_qt[2])
  )

  # Format table with knitr
  print(kable(out, digits = 3, align = "c", col.names = c("Statistic", "Estimate", "CI Lower", "CI Upper")))
  cat("\n")
  cat("Key Decision Rule:\n")
  cat("------------------------------------------------------------\n")
  cat("The null hypothesis should be rejected when Pr(mu_t - mu_c > 0 | Data) > ",
             round(x$opt.gamma$gamma.opt, 5),".\n")
  cat("This threshold accounts for interim looks and represents an increase from", (1 - x$alpha), "\n")
  cat("corresponding to an alpha of", x$alpha*100, "% under a single-look design.\n")
  cat("------------------------------------------------------------\n")

  invisible(out) # Optionally return the table for further processing
}



#' Extract the power of a bcts simulation, with optional adjustment for non-binding futility
#'
#' @param x An object of class bcts
#' @param adjust_for_futility Logical, whether to adjust the power for simulations where non-binding futility was triggered. Default is FALSE.
#' @param conf.level Numeric, the confidence level for the power calculation. Default is 0.95.
#' @param ... Optional arguments
#'
#' @method power bcts
#' @export
#'
#' @importFrom rlang .data
power.bcts <- function(x, adjust_for_futility = FALSE, conf.level = 0.95, ...) {

  if (adjust_for_futility) {
    # Calculate power based on the valid simulations (excluding those with futility-triggering without rejection)

    valid_rows <- which(!x$simresults$fut.trig)

    power <- mc_error_proportion(x = sum(x$simresults$rejectH0.final[valid_rows]),
                                 n = nrow(x$simresults),
                                 level = conf.level)
  } else {
    # Regular power calculation based on rejection of H0 without adjusting for futility
    power <- mc_error_proportion(x = sum(x$simresults$rejectH0.final),
                                 n = nrow(x$simresults),
                                 level = conf.level)
  }

  return(power)
}





#' Evaluate the posterior probability of treatment superiority
#'
#' @param data Trial data
#' @param margin Superiority margin (default: 0)
#' @param gamma Level to declare success (default: 0.975)
#' @param method Method to evaluate superiority (Choose 'mcmc' or 'bayes')
#' @param num_chains The number of parallel chains for the model
#' @param n.adapt The number of iterations for adaptation.
#' @param n.iter 	The number of iterations to monitor
#' @param perc_burnin The percentage of iterations to keep for burn-in
#' @param trt_ref Reference treatment name
#'
#' @description
#' Superiority is established if Pr(mean(treat)-mean(control)>margin) > gamma
#'
#' @references
#' O’Hagan A, Stevens JW, Campbell MJ. Assurance in clinical trial design. Pharmaceut Statist. 2005 Jul;4(3):187–201.
#'
#'
#' @return A dataframe with results from the evaluation
#'
eval_superiority <- function(data,
                             margin, #,
                             gamma, #0.975
                             method,
                             num_chains = 4,
                             n.adapt = 500, n.iter = 1000, perc_burnin = 0.2,
                             trt_ref = "Placebo") {

  treatments <- unique(data$Treatment)
  if (length(treatments) != 2) {
    stop("Superiority can only be assessed for two-arm studies!")
  }

  if (method == "bayes") {
    result <- eval_superiority_bayes(data = data, margin = margin, gamma = gamma,
                                     trt_ref = trt_ref,
                                     num_chains = num_chains, n.adapt = n.adapt,
                                     n.iter = n.iter, perc_burnin = perc_burnin)
  } else if (method == "mcmc") {
    result <- eval_superiority_mcmc(data = data, margin = margin, gamma = gamma, trt_ref = trt_ref)
  } else {
    stop("Unknown method. Please use 'bayes' or 'mcmc'.")
  }

  return(result)
}


#' Bayesian method to evaluate superiority
#'
#' @param data Trial data
#' @param margin Superiority margin (default: 0)
#' @param gamma Level to declare success (default: 0.975)
#' @param trt_ref Reference treatment
#' @param num_chains The number of parallel chains for the model
#' @param n.adapt The number of iterations for adaptation.
#' @param n.iter The number of iterations to monitor
#' @param perc_burnin The percentage of iterations to keep for burn-in
#'
#' @import rjags
#' @importFrom rlang .data
#' @importFrom stats quantile sd qnorm var update
#'
#' @return A data frame with results from the evaluation
#'
eval_superiority_bayes <- function(data, margin, gamma, trt_ref = "Placebo", num_chains, n.adapt, n.iter, perc_burnin) {

  jags.config <- prepare_jags_ppos(dat = data, gamma = gamma, margin = margin, trt_ref = trt_ref)


  # Analyse the interim data
  jags_model <- jags.model(file = textConnection(jags.config$model_string),
                           data = jags.config$jags_data,
                           n.chains = num_chains,
                           n.adapt = n.adapt,
                           quiet = TRUE)

  ##Burnin stage
  update(jags_model, n.iter = ceiling(n.iter*perc_burnin), progress.bar = "none")

  ##Sampling after burnin
  fit_jags_int <- coda.samples(jags_model,
                               variable.names =  jags.config$variable.names,
                               n.iter = floor(n.iter*(1 - perc_burnin)),
                               progress.bar = "none")

  ## Derive posterior distributions at interim
  psample_int <- as.data.frame(do.call(rbind, fit_jags_int))

  mudiff <- psample_int %>% pull(paste0("trteff"))
  mu_t <- psample_int %>% pull(paste0("mu_t"))
  mu_c <- psample_int %>% pull(paste0("mu_c"))
  sigma_t <- psample_int %>% pull(paste0("sigma_t"))
  sigma_c <- psample_int %>% pull(paste0("sigma_c"))
  PPOS <- psample_int %>% pull(paste0("ppos"))

  trt_est <- data.frame("Treatment" = setdiff(data$Treatment, trt_ref),
                        "n_t" = length(!is.na(jags.config$jags_data$Yt)),
                        "n_c" = length(!is.na(jags.config$jags_data$Yc)),
                        "n" = sum(!is.na(data$Y)),
                        "N" = nrow(data),
                        "mu_hat_t" = mean(mu_t),
                        "mu_hat_c" = mean(mu_c),
                        "sigma_hat_t" = mean(sigma_t),
                        "sigma_hat_c" = mean(sigma_c),
                        "est_lower" = quantile(mudiff, 1 - gamma),
                        "est" = mean(mudiff),
                        "est_upper" = quantile(mudiff, gamma),
                        "est_se" = sd(mudiff),
                        "Z_test" = mean(mudiff)/sd(mudiff),
                        "Z_crit" = qnorm(gamma),
                        "I" = 1/var(mudiff),
                        "rejectH0" = quantile(mudiff, (1 - gamma)) > 0,
                        "ppos" = mean(PPOS))

  # Assign the dynamic column names
  names(trt_est)[names(trt_est) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(trt_est)[names(trt_est) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))

  return(trt_est)
}

#' Frequentist method to evaluate superiority
#'
#' @param data Trial data
#' @param margin Superiority margin (default: 0)
#' @param gamma Level to declare success (default: 0.975)
#' @param trt_ref Reference treatment name
#'
#' @return A data frame with the results from the evaluation
#'
#' @importFrom rlang .data
#' @importFrom stats quantile sd qnorm var
eval_superiority_mcmc <- function(data, margin = 0, gamma = 0.975, trt_ref = "Placebo") {

  trt_act <-  setdiff(data$Treatment, trt_ref)
  n <-  sum(!is.na(data$Y))
  N <-  nrow(data)
  Yc <- data %>% dplyr::filter(.data$Treatment == trt_ref & !is.na(.data$Y)) %>% pull("Y")
  Yt <- data %>% dplyr::filter(.data$Treatment == trt_act & !is.na(.data$Y)) %>% pull("Y")

  test_freq <- pow_diff_means(n1 = length(Yt), n2 = length(Yc), N = N,
                              mu1 = mean(Yt), mu2 = mean(Yc), margin = 0,
                              sd1 = sd(Yt), sd2  = sd(Yc),
                              alpha = 1 - gamma,
                              alternative = "superiority")

  trt_est <- data.frame("Treatment" = trt_act,
                        "n_t" = length(Yt),
                        "n_c" = length(Yc),
                        "n" = n,
                        "N" = N,
                        "mu_hat_t" = mean(Yt),
                        "mu_hat_c" = mean(Yc),
                        "sigma_hat_t" = sd(Yt),
                        "sigma_hat_c" = sd(Yc),
                        "est_lower" = test_freq$diff_lower,
                        "est" = test_freq$diff,
                        "est_upper" = test_freq$diff_upper,
                        "est_se" = test_freq$diff_se,
                        "Z_test" = test_freq$Z_test,
                        "Z_crit" = test_freq$Z_crit,
                        "I" = test_freq$I,
                        "rejectH0" = test_freq$rejectH0,
                        "ppos" = test_freq$ppos)

  # Assign the dynamic column names
  names(trt_est)[names(trt_est) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(trt_est)[names(trt_est) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))

  return(trt_est)
}




prepare_jags_ppos <- function(dat, gamma, margin = 0, trt_ref = "Placebo") {

  model_string <- "model {\n"


  model_string <- paste0(model_string, "   for (i in 1:length(Yc)) {\n")
  model_string <- paste0(model_string, "      Yc[i] ~ dnorm(mu_c, prec_c)\n")
  model_string <- paste0(model_string, "   }\n")
  model_string <- paste0(model_string, "   for (i in 1:length(Yt)) {\n")
  model_string <- paste0(model_string, "      Yt[i] ~ dnorm(mu_t, prec_t)\n")
  model_string <- paste0(model_string, "   }\n")
  model_string <- paste0(model_string, "   mu_c ~ dnorm(0, 0.001)\n")
  model_string <- paste0(model_string, "   prec_c ~ dgamma(0.001, 0.001)\n")
  model_string <- paste0(model_string, "   sigma_c <- 1 / sqrt(prec_c)\n\n")
  model_string <- paste0(model_string, "   mu_t ~ dnorm(0, 0.001)\n")
  model_string <- paste0(model_string, "   prec_t ~ dgamma(0.001, 0.001)\n")
  model_string <- paste0(model_string, "   sigma_t <- 1 / sqrt(prec_t)\n\n")

  model_string <- paste0(model_string, "   trteff<- mu_t - mu_c - ", margin, "\n")
  model_string <- paste0(model_string, "   tau <- sqrt(sigma_c**2/length(Yc) + sigma_t**2/length(Yt))\n")
  model_string <- paste0(model_string, "   ppos <- trteff > tau*z_alpha\n\n")
  model_string <- paste0(model_string, "}")

  # Split the data frame by 'trt' and extract the 'Y' values
  Y_list <- split(dat$Y, dat$Treatment)

  # Rename the list elements to Y1, Y2, Y3, etc.
  names(Y_list) <-  ifelse(names(Y_list) == trt_ref, "Yc", "Yt")

  Y_list$z_alpha <- qnorm(gamma)

  variable.names =  c("trteff", "mu_c", "mu_t", "sigma_c", "sigma_t", "ppos")


  return(list(model_string = model_string, jags_data = Y_list,
              variable.names = variable.names))
  }

#' Filter a list of bcts objects based on design criteria
#'
#' This function takes a list of bcts objects and filters them based on user-specified design criteria.
#' The function returns a list of bcts objects that meet the given criteria.
#'
#' @param bcts_list A list of objects of class bcts to be filtered.
#' @param no.looks An optional filter for the number of looks. Only bcts objects with this number of looks will be retained.
#' @param th.fut An optional filter for the futility threshold, checked using approximate equality. Only bcts objects with this threshold will be retained.
#' @param th.eff An optional filter for the efficacy threshold, checked using approximate equality. Only bcts objects with this threshold will be retained.
#' @param n_int1 Optional, the interim sample size in the first stage. Only bcts objects with this interim sample size will be retained.
#' @param n_int2 Optional, the interim sample size in the second stage. Only bcts objects with this interim sample size will be retained.
#' @param n_pln Optional, the planned total sample size. Only bcts objects with this planned sample size will be retained.
#' @param n_max Optional, the maximum allowable sample size. Only bcts objects with this maximum sample size will be retained.
#' @param gamma_values An optional filter for gamma values, checked using approximate equality. Only bcts objects with the specified gamma values will be retained.
#' @param trt_rank An optional named vector specifying treatment ranking. Each bcts object with a matching treatment ranking will be retained.
#' @param power_lower Optional, only retain simulations with power greater than this value (without futility adjustment).
#' @param power_upper Optional, only retain simulations with power less than this value (without futility adjustment).
#' @param select_max_power Logical, if TRUE, only the bcts object with the highest power will be returned.
#' @param tolerance Numeric value for the tolerance when comparing gamma and futility threshold values. Default is 1e-8.
#' @param ... Additional criteria for filtering the bcts objects.
#'
#' @details
#' Filtering by gamma uses approximate equality, meaning that values for certain arguments (e.g., gamma) will match if they are within a small tolerance.
#' Filtering by power is based on the power estimates without any adjustment for futility.
#' If \code{select_max_power} is TRUE, the function will return the bcts object with the highest power
#' after applying the filters.
#'
#' @return A filtered list of bcts objects that meet the design criteria.
#'
#' @export
filter_bcts_by_design <- function(bcts_list,
                                  no.looks = NULL,
                                  th.fut = NULL,
                                  th.eff = NULL,
                                  n_int1 = NULL,
                                  n_int2 = NULL,
                                  n_pln = NULL,
                                  n_max = NULL,
                                  gamma_values = NULL,
                                  trt_rank = NULL,
                                  power_lower = NULL,
                                  power_upper = NULL,
                                  select_max_power = FALSE,
                                  tolerance = 1e-8, ...) {
  # Filter by number of looks
  if (!is.null(no.looks)) {
    bcts_list <- Filter(function(x) x$no.looks == no.looks, bcts_list)
  }

  # Filter by futility threshold using approximate equality
  if (!is.null(th.fut)) {
    bcts_list <- Filter(function(x) {
      abs(x$th.fut - th.fut) < tolerance
    }, bcts_list)
  }

  # Filter by futility threshold using approximate equality
  if (!is.null(th.eff)) {
    bcts_list <- Filter(function(x) {
      abs(x$th.eff - th.eff) < tolerance
    }, bcts_list)
  }

  # Filter by interim sample size for the first stage
  if (!is.null(n_int1)) {
    bcts_list <- Filter(function(x) {
      !is.null(x$n$Int1) && abs(x$n$Int1 - n_int1) < tolerance
    }, bcts_list)
  }

  # Filter by interim sample size for the second stage
  if (!is.null(n_int2)) {
    bcts_list <- Filter(function(x) {
      !is.null(x$n$Int2) && abs(x$n$Int2 - n_int2) < tolerance
    }, bcts_list)
  }

  # Filter by planned sample size
  if (!is.null(n_pln)) {
    bcts_list <- Filter(function(x) {
      !is.null(x$n$Planned) && abs(x$n$Planned - n_pln) < tolerance
    }, bcts_list)
  }

  # Filter by maximum size
  if (!is.null(n_max)) {
    bcts_list <- Filter(function(x) {
      !is.null(x$n$Maximum) && abs(x$n$Maximum - n_max) < tolerance
    }, bcts_list)
  }

  # Filter by gamma values using approximate equality
  if (!is.null(gamma_values)) {
    bcts_list <- Filter(function(x) {
      any(abs(x$gamma - gamma_values) < tolerance)
    }, bcts_list)
  }

  # Filter by treatment rank
  if (!is.null(trt_rank)) {
    bcts_list <- Filter(function(x) {
      if (is.null(x$trt_rank)) {
        return(FALSE)  # Skip if trt_rank is NULL in the object
      }
      all(names(trt_rank) %in% names(x$trt_rank) & trt_rank == x$trt_rank)
    }, bcts_list)
  }

  # Filter by power (lower and upper bounds)
  if (!is.null(power_lower)) {
    bcts_list <- Filter(function(x) {
      xpow <- power(x, adjust_for_futility = FALSE)
      return(xpow$est > power_lower)
      } , bcts_list)
  }

  if (!is.null(power_upper)) {
    bcts_list <- Filter(function(x) {
      xpow <- power(x, adjust_for_futility = FALSE)
      return(xpow$est < power_upper)
      }, bcts_list)
  }

  # Check if the list is empty after filtering
  if (length(bcts_list) == 0) {
    return(NULL)  # Return NULL if no objects match the criteria
  }

  # If select_max_power is TRUE, return the bcts object with the highest power
  if (select_max_power && length(bcts_list) > 0) {
    max_power_object <- bcts_list[[which.max(sapply(bcts_list, function(x) power(x, adjust_for_futility = FALSE)$est))]]
    return(max_power_object)
  }


  return(bcts_list)
}


