#' Bayesian Clinical Trial Simulation (BCTS)
#'
#' @description
#' This function performs Bayesian clinical trial simulations to evaluate adaptive trial designs.
#' It allows for up to two interim analyses, one for dose selection and one for sample size re-estimation.
#' The null hypothesis is rejected if the posterior probability \eqn{Pr(\mu_t - \mu_c > 0 | \text{Data}) > \gamma} exceeds a calibrated threshold.
#'
#' @param n_dose_sel Sample size at the first interim analysis for dose selection.
#' @param n_ss_reest Optional sample size at the second interim analysis for sample size re-estimation.
#'                   If not provided, only one interim analysis is performed.
#' @param n_pln Planned sample size for the trial.
#' @param n_max Maximum allowable sample size for the trial.
#' @param mu Named vector of expected effect sizes for each treatment group (including the control group).
#' @param sigma Named vector of expected standard deviations for each treatment group.
#' @param trt_ref Character string specifying the control treatment (reference group)
#' @param trt_rank Named vector of preference rankings for each treatment group.
#'                 Lower values indicate higher preference, used when multiple treatments exceed the efficacy threshold (`th.eff`).
#' @param alpha Desired one-sided Type-I error rate. Default is `0.025`.
#'              This represents the maximum allowable probability of rejecting the null hypothesis when it is true.
#' @param alpha_tolerance Tolerance for the Type-I error rate. Default is `0.001`.
#'                        This specifies the acceptable deviation from the target `alpha` during the calibration process.
#'                        Smaller values indicate stricter adherence to the target, while larger values allow more flexibility.
#' @param th.fut Predictive power threshold for declaring futility. Default is `0.2`.
#' @param th.eff Predictive power threshold for declaring efficacy. Default is `0.9`.
#' @param th.prom  Predictive power threshold that triggers sample size re-estimation. Default is `0.5`.
#' @param method Method for power estimation. Options are `"bayes"` for a fully Bayesian approach or `"mcmc"` for a frequentist approximation. Default is `"mcmc"`.
#' @param nsim Number of trials to simulate (default: 1000)
#' @param num_chains Number of MCMC chains (default: 4)
#' @param n.iter Total number of MCMC iterations. Default is `5000`.
#' @param n.adapt Number of adaptation iterations for MCMC. Default is `500`.
#' @param perc_burnin Proportion of MCMC iterations used for burn-in. Default is `0.2`.
#' @param progress.bar Type of progress bar displayed during simulations. Options are `"text"`, `"gui"`, and `"none"`.
#'
#' @details
#' The function supports two designs:
#' - **One interim analysis**: Used for dose selection only.
#' - **Two interim analyses**: The first for dose selection and the second for sample size re-estimation.
#'
#' The function begins by calibrating the \eqn{\gamma} threshold to achieve the target Type-I error rate (\eqn{\alpha}).
#' After calibration, the Type-I error and power are assessed under the calibrated threshold.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @return
#' An object of class `"bcts_results"` containing the following components:
#' - **sim_type1**: Results from Type-I error simulations.
#' - **sim_power**: Results from power simulations.
#' - **design**: Details of the trial design, including `alpha`.
#' - **estimation**: A list detailing the estimation and simulation settings:
#'     - `nsim`: Number of simulated trials.
#'     - `method`: Method used for power estimation (`"bayes"` or `"mcmc"`).
#'     - `num_chains`: Number of MCMC chains.
#'     - `n.iter`: Total number of MCMC iterations.
#'     - `n.adapt`: Number of adaptation iterations for MCMC.
#'     - `perc_burnin`: Proportion of MCMC iterations used for burn-in.
#' - **result**: A list with the calibrated results:
#'     - `gamma_threshold`: Calibrated threshold for rejecting the null hypothesis.
#'     - `critical_z_value`: Critical value to evaluate the null hypothesis (\eqn{z = \Phi^{-1}(\gamma)}).
#'     - `critical_p_value`: Critical p-value corresponding to \eqn{z.crit}, representing the smallest significance level at which the null hypothesis should be rejected.
#' - **gamma.table**: A data frame summarizing the gamma calibration process,
#'   including the tested gamma values and their corresponding Type-I error rates.
#'
#' @examples
#' \dontrun{
#' bcts(n_dose_sel = 60, n_ss_reest = 175, n_pln = 200, n_max = 250,
#'      mu = c("Placebo" = 0, "Dose 1" = 0.4, "Dose 2" = 0.5),
#'      sigma = c("Placebo" = 1, "Dose 1" = 1, "Dose 2" = 1),
#'      trt_ref = "Placebo", method = "mcmc", alpha = 0.01,
#'      th.fut = 0.15, th.eff = 0.9, th.prom = 0.5, nsim = 10000)
#' }
#'
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

  # Extract treatment names
  trt_names <- extract_treatment_names(mu = mu, sigma = sigma)

  if (is.null(n_ss_reest)) {
    # Set no.looks
    no.looks <- 1

    # Calibrate gamma for a one-interim design
    message("Using bcts_one_interim.")
    opt.gamma <- calibrate_gamma(
      bcts_fun = function(x, nsim) {
        bcts_one_interim(n_int = n_dose_sel, n_pln = n_pln, n_max = n_max,
                         mu = mu * 0, sigma = sigma, trt_ref = trt_ref,
                         gamma = x,
                         trt_rank = trt_rank, th.fut = th.fut, th.eff = th.eff,
                         th.prom = th.prom, method = method,
                         nsim = nsim,
                         num_chains = num_chains, n.iter = n.iter,
                         n.adapt = n.adapt, perc_burnin = perc_burnin,
                         progress.bar = progress.bar)
      }, type1 = alpha, type1.tolerance = alpha_tolerance, nsim = nsim)

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
    # Set no.looks
    no.looks <- 2

    # Calibrate gamma for the two-interim version
    message("Using bcts_two_interim.")
    opt.gamma <- calibrate_gamma(
      bcts_fun = function(x, nsim) {
        bcts_two_interim(n_int1 = n_dose_sel, n_int2 = n_ss_reest,
                         n_pln = n_pln, n_max = n_max,
                         mu = mu * 0, sigma = sigma, trt_ref = trt_ref,
                         gamma = x, trt_rank = trt_rank, th.fut = th.fut,
                         th.eff = th.eff, th.prom = th.prom, method = method,
                         nsim = nsim,
                         num_chains = num_chains, n.iter = n.iter,
                         n.adapt = n.adapt, perc_burnin = perc_burnin,
                         progress.bar = progress.bar)
      }, type1 = alpha, type1.tolerance = alpha_tolerance,nsim = nsim)

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

  # Derive summary metrics
  z_crit <- qnorm(opt.gamma$gamma.opt)

  out <- list(sim_type1 = opt.gamma$sim.opt,
              sim_power = sim_power,
              design = list(no.looks = no.looks,
                            alpha = alpha,
                            trt_ref = trt_ref,
                            trt_active = setdiff(trt_names, trt_ref),
                            trt_rank = trt_rank,
                            mu = mu,
                            sigma = sigma,
                            th.fut = th.fut,
                            th.eff = th.eff,
                            th.prom = th.prom),
              estimation = list(nsim = nsim,
                                method = method,
                                num_chains = num_chains,
                                n.iter = n.iter,
                                n.adapt = n.adapt,
                                perc_burnin = perc_burnin),
              result = list(
                type1_error = power(opt.gamma$sim.opt, adjust_for_futility = FALSE),                           # Clearly indicates Type-I error
                power_no_futility = power(sim_power, adjust_for_futility = FALSE),         # Describes power without futility adjustment
                power_with_futility = power(sim_power, adjust_for_futility = TRUE),          # Describes power with futility adjustment
                gamma_threshold = opt.gamma$gamma.opt,         # Indicates the gamma threshold
                critical_z_value = z_crit,                     # Clearly indicates the critical z-value
                critical_p_value = pnorm(z_crit, lower.tail = FALSE) # Describes the critical p-value
              ),
              gamma.table = opt.gamma$table)
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
#' @importFrom purrr map_chr
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

  ## Header
  cat("Sample Size Calculation for an Adaptive Trial\n\n")
  cat("Design Summary:\n")
  cat("- Target Type-I Error:", sprintf("%.2f%%", x$design$alpha * 100), "\n")
  cat("- Number of Looks:", x$design$no.looks, "\n\n")

  ## PPOS Thresholds
  cat("PPOS Thresholds:\n")
  cat("- Futility:", sprintf("%.2f%%", x$design$th.fut * 100), "\n")
  cat("- Efficacy:", sprintf("%.2f%%", x$design$th.eff * 100), "\n")
  cat("- Sample Size Re-estimation:", sprintf("%.2f%%", x$design$th.prom * 100), "\n\n")

  ## Effect Sizes and SD
  cat("Expected Effect Sizes and Standard Deviations:\n")
  effect_sizes <- union(x$design$trt_ref, x$design$trt_active) %>%
    map_chr(function(treatment) {
      effect_size <- ifelse(is.null(x$sim_power$mu[[treatment]]), "N/A", x$sim_power$mu[[treatment]])
      sd <- ifelse(is.null(x$sim_power$sigma[[treatment]]), "N/A", x$sim_power$sigma[[treatment]])
      is_reference <- ifelse(treatment == x$design$trt_ref, " (Reference)", "")
      paste("-", treatment, ": Effect Size =", effect_size, ", SD =", sd, is_reference)
    })
  cat(paste(effect_sizes, collapse = "\n"), "\n\n")

  ## Sample Size Details
  cat("Sample Size Details:\n")
  if (x$design$no.looks == 2) {
    cat("- Interim 1:", x$sim_power$n$Int1, "patients with outcome data\n")
    cat("- Interim 2:", x$sim_power$n$Int2, "patients with outcome data\n")
  } else if (x$design$no.looks == 2) {
    cat("- Interim:", x$sim_power$n$Int1, "patients with outcome data\n")
  }
  cat("- Planned Sample Size:", x$sim_power$n$Planned, "patients\n")
  cat("- Maximum Sample Size:", x$sim_power$n$Maximum, "patients\n\n")

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

  if (x$estimation$method == "mcmc") {
    # Provide critical p-value and z-value for mcmc method
    cat("The null hypothesis should be rejected when the critical p-value is below",
        sprintf("%.5f", x$result$critical_p_value), "or equivalently when the critical z-value exceeds",
        sprintf("%.3f", x$result$critical_z_value), ".\n")
    cat("This threshold accounts for interim looks and represents an increase from the single-look z-value of",
        sprintf("%.3f", qnorm(1 - x$design$alpha)), ".\n")
  } else {
    # Provide gamma threshold for bayes method
    cat("The null hypothesis should be rejected when Pr(mu_t - mu_c > 0 | Data) >",
        sprintf("%.5f", x$result$gamma), ".\n")
    cat("This threshold accounts for interim looks and represents an increase from",
        sprintf("%.5f", 1 - x$design$alpha), ".\n")
    cat("corresponding to an alpha of", sprintf("%.2f%%", x$design$alpha * 100),
        "under a single-look design.\n")
  }
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
                              mu1 = mean(Yt), mu2 = mean(Yc), margin = margin,
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

#' Generate a PDF Report for a bcts_results Object
#'
#' @param x An object of class `bcts_results` containing the results to be included in the report.
#' @param output_file A character string specifying the output PDF file name (default: "sample_size_report.pdf").
#' @param output_dir A character string specifying the directory where the PDF file will be saved (default: current working directory).
#' @param template_path A character string specifying the path to a custom Quarto template file. If not provided, the default template in the package is used.
#' @return The path to the generated PDF report.
#' @export
generate_report <- function(x, output_file = "sample_size_report.pdf", output_dir = getwd(), template_path = NULL) {
  # Check if Quarto is available
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("The 'quarto' package is required to generate the report. Please install it using install.packages('quarto').")
  }

  # Validate bcts_results
  if (is.null(x)) {
    stop("The 'x' object cannot be NULL.")
  }

  # Determine the template path
  if (is.null(template_path)) {
    template_path <- system.file("rmarkdown/templates/report.qmd", package = "bcts")
    if (template_path == "") {
      stop("Default Quarto template not found. Ensure that the 'report.qmd' file is included in your package under 'inst/rmarkdown/templates/'.")
    }
  } else {
    if (!file.exists(template_path)) {
      stop("The specified custom template file does not exist: ", template_path)
    }
  }

  # Generate the report in the current working directory
  temp_output_file <- output_file  # File name only
  quarto::quarto_render(
    input = template_path,
    output_format = "pdf",
    output_file = temp_output_file,
    execute_params = list(sim = x)
  )

  # Move the file to the specified output directory
  final_output_path <- file.path(output_dir, output_file)
  file.rename(temp_output_file, final_output_path)

  # Return the final path to the report
  return(final_output_path)
}

