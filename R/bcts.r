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
#' @importFrom stats pnorm
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


#' Bayesian clinical trial simulation (BCTS)
#'
#' @param n_int Interim sample size
#' @param n_pln Planned sample size
#' @param n_max Maximum sample size
#' @param mu Named vector with expected effect sizes
#' @param sigma Named vector with expected standard deviations
#' @param trt_ref Character denoting the control treatment
#' @param trt_rank Named vector with preference ranking for each treatment (e.g., lower doses are preferred)
#' @param prioritize_low_rank  If multiple treatments have a posterior predictive power > 'th.eff', the treatment with lowest rank
#' will be selected at interim
#' @param gamma Level to declare success. For a fixed design, gamma is typically chosen as 0.975 for a one-sided type-I error rate of 2.5%. However, an increase is usually needed because of dose selection.
#' @param th.fut Futility threshold for the predictive power calculated conditioned on the interim data
#' @param th.eff Efficacy threshold for the predictive power calculated conditioned on the interim data
#' @param th.prom Predictive power threshold that would trigger sample size increase
#' @param method Method to estimate the power. Choose "bayes" for a fully Bayesian approach or "mcmc" for a frequentist approximation.
#' @param nsim Number of trials to simulate (default: 1000)
#' @param num_chains Number of MCMC chains
#' @param n.iter Number of MCMC iterations for estimation
#' @param n.adapt Number of MCMC iterations for adaptation
#' @param perc_burnin How many n.iter should be reserved for burnin?
#' @param progress.bar Type of progress bar. Possible values are "text", "gui", and "none"
#'
#' @description
#' The null hypothesis can be rejected when the posterior probability that $\\mu_t$ - $\\mu_c$ exceeds a high probability threshold,
#' i.e. $Pr(\\mu_t - \\mu_c > 0|D) > $\\gamma$)$.
#'
#'
#' @return An object of class "bcts"
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr slice_head mutate select filter add_row pull
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom binom binom.confint
#' @importFrom rlang .data
#'
bcts <- function(n_int, n_pln, n_max, mu, sigma, trt_ref, trt_rank,
                 prioritize_low_rank = TRUE, gamma = 0.975,
                 th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                 method = "mcmc", nsim = 1000, num_chains = 4,
                 n.iter = 5000, n.adapt = 500, perc_burnin = 0.2,
                 progress.bar = "text") {

  # Sample the random seeds
  seeds <- seq(nsim)

  n_arms <- length(mu)

  # Extract the treatment names
  trt_names <- names(mu)
  names(mu) <- names(sigma) <- trt_names

  # Active treatment names
  trt_active <- setdiff(trt_names, trt_ref)

  simresults <- data.frame(sim = seq(nsim),
                           seed = NA,
                           fut.trig = FALSE,
                           inc.ss = FALSE,
                           sel.dose = NA,
                           sel.dose.ppos = NA,
                           n.final = NA,
                           est.final = NA,
                           est_lower.final = NA,
                           est_se.final = NA,
                           est_upper.final = NA,
                           rejectH0.final = NA)

  n.interim <- n.final <- data.frame(sim = numeric(), Treatment = character(), n = integer(), stringsAsFactors = FALSE)


  if (progress.bar == "text") {
    pb <- utils::txtProgressBar(min = 0, max = nsim, initial = 0)
  }

  for (i in 1:nsim) {
    simresults$seed[i] <- seeds[i]
    set.seed(seeds[i])

    # Simulate trial data at interim stage
    dat_int <- sim_rct_normal(n = n_int,
                              mean = mu,
                              sd = sigma,
                              trtnames = trt_names)

    # Summarize the data
    n_int_summary_data <- dat_int %>%
      dplyr::group_by(.data$Treatment) %>%
      dplyr::summarize(n = dplyr::n(), .groups = 'drop') %>%
      dplyr::mutate(sim = i)
    n.interim <- dplyr::bind_rows(n.interim , n_int_summary_data)

    ri <- dose_selection(dat_int = dat_int, n_pln = n_pln, n_max = n_max,
                         trt_ref = trt_ref,
                              trt_active = trt_active, trt_rank = trt_rank,
                              prioritize_low_rank = prioritize_low_rank,
                              gamma = gamma, th.fut = th.fut, th.eff = th.eff,
                              th.prom = th.prom,
                              method = method, num_chains = num_chains,
                              n.iter = n.iter, n.adapt = n.adapt, perc_burnin = perc_burnin)

    simresults$fut.trig[i] <- ri$result$fut.trig
    simresults$inc.ss[i] <- ri$result$inc.ss
    simresults$n.final[i] <- ri$result$n.final
    simresults$sel.dose[i] <- ri$result$sel.dose
    simresults$sel.dose.ppos[i] <- ri$result$sel.dose.ppos

    # Efficacy analysis in the final dataset
    dat_xtr <- sim_rct_normal(n = simresults$n.final[i] - nrow(dat_int),
                              mean = mu[c(trt_ref, simresults$sel.dose[i])],
                              sd = sigma[c(trt_ref, simresults$sel.dose[i])],
                              trtnames = c(trt_ref, simresults$sel.dose[i]))

    # Summarize the data
    n_fin_summary_data <- rbind(dat_int, dat_xtr) %>%
      dplyr::group_by(.data$Treatment) %>%
      dplyr::summarize(n = dplyr::n(), .groups = 'drop') %>%
      dplyr::mutate(sim = i)
    n.final <- dplyr::bind_rows(n.final , n_fin_summary_data)

    dat_fin <- rbind(dat_int %>% dplyr::filter(.data$Treatment %in% c(trt_ref, simresults$sel.dose[i])), dat_xtr)
    dat_fin$Treatment <- droplevels(dat_fin$Treatment)

    result_fin <- eval_superiority(data = dat_fin, margin = 0, gamma = gamma,
                                   trt_ref = trt_ref, method = method,
                                   num_chains = num_chains, n.adapt = n.adapt,
                                   n.iter = n.iter, perc_burnin = perc_burnin)

    simresults$est.final[i] <- result_fin$est
    simresults$est_lower.final[i] <-  result_fin[paste0("est_", sub("0\\.", "", 1 - gamma))]
    simresults$est_upper.final[i] <-  result_fin[paste0("est_", sub("0\\.", "", gamma))]
    simresults$est_se.final[i] <- result_fin$est_se
    simresults$rejectH0.final[i] <- result_fin$rejectH0

    if (progress.bar == "text") {
      utils::setTxtProgressBar(pb, i)
    }
  }

  if (progress.bar == "text") {
    close(pb)
  }

  names(simresults)[names(simresults) == "est_lower.final"] <- paste0("est_", sub("0\\.", "", 1 - gamma), ".final")

  out <- list(trt_ref = trt_ref,
              mu = mu,
              sigma = sigma,
              gamma = gamma,
              nsim = nsim,
              method = method,
              n.interim = n.interim,
              n.final = n.final,
              simresults = simresults
  )

  class(out) <- "bcts"
  return(out)
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

#' Extract the power
#'
#' @param \dots Optional arguments
#'
#' @details This is a generic function.
#'
#' @export power
power <- function(...) {
  UseMethod("power")
}


#' Plot selected dose
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotSelDose
plotSelDose <- function(...) {
  UseMethod("plotSelDose")
}

#' Plot the final sample size
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotFinalSampleSize
plotFinalSampleSize  <- function(...) {
  UseMethod("plotFinalSampleSize")
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

  n_fin_by_arm <- x$n.final %>%
    dplyr::group_by(.data$Treatment) %>%
    dplyr::summarize(est = mean(.data$n),
                     cil = quantile(.data$n, (1 - conf.level)/2),
                     ciu = quantile(.data$n, 1 - (1 - conf.level)/2)) %>%
    dplyr::mutate(statistic = paste0("N final (", .data$Treatment, ")")) %>%
    dplyr::select(-"Treatment")

  n_select_arm <- x$simresults %>%
    dplyr::group_by(.data$sel.dose) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::mutate(statistic = paste0("Pr(select ", .data$sel.dose, ")"),
                  est = .data$n/nrow(x$simresults),
                  cil = binom.confint(x = .data$n,
                                        n = nrow(x$simresults),
                                        methods = "exact", conf.level = conf.level)$lower,
                    ciu = binom.confint(x = .data$n,
                                        n = nrow(x$simresults),
                                        methods = "exact", conf.level = conf.level)$upper) %>%
    dplyr::select(-c("sel.dose", "n"))

  power <- mc_error_proportion(x = sum(x$simresults$rejectH0.final),
                               n = nrow(x$simresults),
                               level = conf.level)

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
  out <- out %>% add_row(data.frame(statistic = "Power",
                                      est = power$est,
                                      cil = power$lower,
                                      ciu = power$upper))
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
    out <- out %>% add_row(n_fin_by_arm)
    out <- out %>% add_row(n_select_arm)
    print(out)
}



#' Plot simulation results
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom rlang .data
plot.bcts <- function(x, ...) {
  conf.level <- 0.95

  power <- mc_error_proportion(x = sum(x$simresults$rejectH0.final),
                               n = nrow(x$simresults),
                               level = conf.level)

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

  out <- out %>% add_row(data.frame(statistic = "Power",
                                    est = power$est,
                                    cil = power$lower,
                                    ciu = power$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(fut.trig)",
                                    est = fut.trig$est,
                                    cil = fut.trig$lower,
                                    ciu = fut.trig$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(inc.ss)",
                                    est = inc.ss$est,
                                    cil = inc.ss$lower,
                                    ciu = inc.ss$upper))

  # Convert to percentages
  out$est <- out$est * 100
  out$cil <- out$cil * 100
  out$ciu <- out$ciu * 100

  ggplot(out, aes(x = .data$statistic, y = .data$est)) +
    #geom_bar(stat = "identity", alpha = 0.7) +
    geom_point() +
    geom_errorbar(aes(ymin = .data$cil, ymax = .data$ciu), width = 0.2) +
    facet_wrap(~statistic, scales = "free") +
    labs(x = "Statistic", y = "Estimate (95% CI)") +
    xlab("") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) #+ scale_fill_brewer()

}

#' Plot posterior distribution of the mean outcome
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
posterior.bcts <- function(x, ...) {
  ggplot(x$simresults, aes(x = .data$est.final)) +
    geom_density() +
    xlab("Mean outcome selected dose")
}

#' Extract the power of a bcts simulation
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @importFrom rlang .data
power.bcts <- function(x, ...) {
  conf.level <- 0.95

  power <- mc_error_proportion(x = sum(x$simresults$rejectH0.final),
                               n = nrow(x$simresults),
                               level = conf.level)
  return(power)
}


#' Plot selected dose
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotSelDose.bcts <- function(x, ...) {

  # Calculate the percentage for each bar
  nDose <- x$simresults %>%
    dplyr::count(.data$sel.dose) %>%
    mutate(percentage = .data$n / sum(.data$n) * 100)

  # Plot the bar chart with percentage labels
  ggplot(nDose, aes(x = .data$sel.dose, y = .data$n, fill = .data$sel.dose)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(.data$percentage, 1), "%")), vjust = -0.5) +
    labs(x = "Selected Dose", y = "Count") +
    scale_fill_brewer()  +
    theme(legend.position = "none") +
    ylim(0,  nrow(x$simresults))

}

#' Plot final sample size
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotFinalSampleSize.bcts <- function(x, ...) {

  # Join n.final with selected treatment
  #dat <- x$n.final %>% merge(x$simresults %>% select("sim", "sel.dose"))

  #ggplot(dat, aes(x=n, group = Treatment, fill = Treatment)) +
  #  geom_histogram(col = "black") +
  #  facet_wrap(~sel.dose) + scale_fill_brewer()

  # Calculate the mean final sample size
  mean_n_final <- mean(x$simresults$n.final)

  ggplot(x$simresults, aes(x = .data$n.final)) +
    geom_histogram() +
    xlab("Final sample size") +
    ylab("Count") +
    annotate("text", x = Inf, y = Inf, label = paste("Mean final sample size:", round(mean_n_final, 0)),
             hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +
    ylim(0,  nrow(x$simresults))
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
#' @import stats
#' @import rjags
#' @importFrom rlang .data
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

#' Apply dose selection at interim
#'
#' @param dat_int Observed data at the interim stage
#' @param n_pln Planned sample size
#' @param n_max Maximum sample size
#' @param trt_ref Character denoting the control treatment
#' @param trt_active Character vector denoting the active treatment names
#' @param trt_rank Named vector with preference ranking for each treatment (e.g., lower doses are preferred)
#' @param prioritize_low_rank  If multiple treatments have a posterior predictive power > 'th.eff', the treatment with lowest rank
#' will be selected at interim
#' @param gamma Level to declare success. For a fixed design, gamma is typically chosen as 0.975 for a one-sided type-I error rate of 2.5%. However, an increase is usually needed because of dose selection.
#' @param th.fut Futility threshold for the predictive power calculated conditioned on the interim data
#' @param th.eff Efficacy threshold for the predictive power calculated conditioned on the interim data
#' @param th.prom Predictive power threshold that would trigger sample size increase
#' @param method Method to estimate the power. Choose "bayes" for a fully Bayesian approach or "mcmc" for a frequentist approximation.
#' @param num_chains Number of MCMC chains
#' @param n.iter Number of MCMC iterations for estimation
#' @param n.adapt Number of MCMC iterations for adaptation
#' @param perc_burnin How many n.iter should be reserved for burnin?
#'
#' @return A list with PPoS and a dataframe containing information on the selected dose and final sample size
#'
#' @export
#'
#' @importFrom rlang .data
dose_selection <- function(dat_int, n_pln, n_max, trt_ref, trt_active,
                           trt_rank,
                                prioritize_low_rank = TRUE, gamma = 0.975,
                                th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                                method = "mcmc", num_chains = 4,
                                n.iter = 5000, n.adapt = 500, perc_burnin = 0.2) {

  # Evaluate PPOS for each treatment at the planned sample size
  ppos <- setNames(rep(NA, length(trt_active)), trt_active)
  result_final <- data.frame(fut.trig = NA, inc.ss = NA, n.final = NA, sel.dose = "", sel.dose.ppos = NA)

  for (trt_act_i in trt_active) {

    dat_xtr <- sim_rct_normal(n = n_pln - nrow(dat_int),
                              mean = rep(NA,2),
                              sd = rep(NA,2),
                              trtnames = c(trt_ref, trt_act_i))
    dat_pln <- rbind(dat_int %>% dplyr::filter(.data$Treatment %in% c(trt_ref, trt_act_i)), dat_xtr)
    dat_pln$Treatment <- droplevels(dat_pln$Treatment)

    result_pln <- eval_superiority(dat_pln, margin = 0, gamma = gamma,
                                   method = method, trt_ref = trt_ref,
                                   num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)
    ppos[trt_act_i] <-  result_pln %>% pull("ppos")
  }

  if (max(ppos) < th.fut) {
    result_final$fut.trig <- TRUE
    result_final$inc.ss <- FALSE
    result_final$n.final <- n_pln
    result_final$sel.dose <- names(ppos)[which.max(ppos)]
    result_final$sel.dose.ppos <- max(ppos)
  } else if (max(ppos) >= th.prom & max(ppos) < th.eff) {
    result_final$fut.trig <- FALSE
    result_final$inc.ss <- TRUE

    # Evaluate predictive power for each sample size
    ppos_by_n <- data.frame(n_per_arm_new = numeric(),
                            Treatment = character(),
                            n = numeric(),
                            N = numeric(),
                            ppos = numeric())

    for (n_pln_new in (n_pln + 1):n_max) {
      ppos_eval <- setNames(rep(NA, length(trt_active)), trt_active)
      for (trt_act_i in trt_active) {
        dat_xtr <- sim_rct_normal(n = n_pln_new - nrow(dat_int),
                                  mean = rep(NA, 2),
                                  sd = rep(NA, 2),
                                  trtnames = c(trt_ref, trt_act_i))
        dat_pln <- rbind(dat_int %>% dplyr::filter(.data$Treatment %in% c(trt_ref, trt_act_i)), dat_xtr)
        dat_pln$Treatment <- droplevels(dat_pln$Treatment)

        result_pln <- eval_superiority(dat_pln, margin = 0, gamma = gamma,
                                       method = "mcmc", trt_ref = trt_ref)

        ppos_eval[trt_act_i] <- result_pln %>% pull("ppos")
      }
      result_final$n.final <- n_pln_new
      result_final$sel.dose <- names(ppos_eval)[which.max(ppos_eval)]
      result_final$sel.dose.ppos <- max(ppos_eval)
      if (max(ppos_eval) >= th.eff) {
        # Required sample size increase reached
        break;
      }
    }
  } else {
    result_final$fut.trig <- FALSE
    result_final$inc.ss <- FALSE
    result_final$n.final <- n_pln

    # Filter treatments that have ppos values above th.eff
    filtered_treatments <- ppos[ppos >= th.eff]

    if (length(filtered_treatments) >= 2 & prioritize_low_rank) {
      # Get the treatment name with the lowest rank among the filtered treatments
      treatment_with_lowest_rank <- names(which.min(trt_rank[names(filtered_treatments)]))
      result_final$sel.dose <- treatment_with_lowest_rank
      result_final$sel.dose.ppos <- ppos[treatment_with_lowest_rank]
    } else {
      result_final$sel.dose <- names(ppos)[which.max(ppos)]
      result_final$sel.dose.ppos <- max(ppos)
    }
  }


  return(list(PPoS = ppos, result = result_final))
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



