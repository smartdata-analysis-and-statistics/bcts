search_gamma <- function() {

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
#' @param futility_threshold An optional filter for the futility threshold (e.g., 0.15, 0.20). Only bcts objects with this threshold will be retained.
#' @param gamma_values An optional filter for gamma values. Only bcts objects with the specified gamma values will be retained.
#' @param ... Additional criteria for filtering the bcts objects.
#'
#' @return A filtered list of bcts objects that meet the design criteria.
#'
#' @export
filter_bcts_by_design <- function(bcts_list, no.looks = NULL, futility_threshold = NULL, gamma_values = NULL, ...) {
  # Filter by number of looks
  if (!is.null(no.looks)) {
    bcts_list <- Filter(function(x) x$no.looks == no.looks, bcts_list)
  }

  # Filter by futility threshold
  if (!is.null(futility_threshold)) {
    bcts_list <- Filter(function(x) x$th.fut == futility_threshold, bcts_list)
  }

  # Filter by gamma values
  if (!is.null(gamma_values)) {
    bcts_list <- Filter(function(x) x$gamma %in% gamma_values, bcts_list)
  }

  return(bcts_list)
}


