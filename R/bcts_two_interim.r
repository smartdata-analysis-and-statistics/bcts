#' Bayesian Clinical Trial Simulation with Two Interim Analyses (BCTS)
#'
#' @param n_int First interim sample size (for dose selection)
#' @param n_int2 Second interim sample size (for sample size re-estimation)
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
#' This function simulates Bayesian clinical trials with two interim analyses:
#' the first for dose selection and the second for sample size re-estimation.
#' The null hypothesis can be rejected when the posterior probability that
#' $\\mu_t$ - $\\mu_c$ exceeds a high probability threshold, i.e.,
#' $Pr(\\mu_t - \\mu_c > 0|D) > \\gamma$.
#'
#'
#' @export
#' @return An object of class "bcts"
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr slice_head mutate select filter add_row pull
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom binom binom.confint
#' @importFrom rlang .data
#'
bcts_two_interim <- function(n_int, n_int2, n_pln, n_max, mu, sigma, trt_ref, trt_rank,
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

  # Keep track of futility triggering (each column refers to an interim analysis)
  fut.trig <- data.frame(matrix(NA, nrow = nsim, ncol = 2))

  # Keep track of treatment benefit
  int1.ben <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = length(trt_active))), trt_active) #interim 1
  int2.ben <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = length(trt_active))), trt_active) #interim 2

  # Keep track of PPoS
  int1.PPoS <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = length(trt_active))), trt_active) #interim 1
  int2.PPoS <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = length(trt_active))), trt_active) #interim 2

  # Keep track of selected dose at first interim analysis
  sel.dose <- rep("", nsim)

  if (progress.bar == "text") {
    pb <- utils::txtProgressBar(min = 0, max = nsim, initial = 0)
  }

  for (i in 1:nsim) {
    simresults$seed[i] <- seeds[i]
    set.seed(seeds[i])

    # Simulate trial data at the first interim
    dat_int <- sim_rct_normal(n = n_int,
                              mean = mu,
                              sd = sigma,
                              trtnames = trt_names)

    # Summarize the data
    n_int1_summary_data <- dat_int %>%
      dplyr::group_by(.data$Treatment) %>%
      dplyr::summarize(n = dplyr::n(), .groups = 'drop') %>%
      dplyr::mutate(sim = i)
    n.int1 <- dplyr::bind_rows(n.interim , n_int1_summary_data)

    ri <- dose_selection(dat_int = dat_int, n_pln = n_pln, n_max = n_max,
                         trt_ref = trt_ref,
                         trt_active = trt_active, trt_rank = trt_rank,
                         prioritize_low_rank = prioritize_low_rank,
                         gamma = gamma, th.fut = th.fut, th.eff = th.eff,
                         th.prom = th.prom,
                         method = method, num_chains = num_chains,
                         n.iter = n.iter, n.adapt = n.adapt, perc_burnin = perc_burnin)

    fut.trig[i,1] <- ri$result$fut.trig
    sel.dose[i] <- ri$result$sel.dose

    # Efficacy at interim 1
    int1.ben[i,] <- ri$benefit

    # PPoS at interim 1
    int1.PPoS[i, ] <- ri$PPoS

    # Continue with interim analysis 2
    dat_xtr <- sim_rct_normal(n = n_int2 - nrow(dat_int),
                              mean = mu[c(trt_ref, sel.dose[i])],
                              sd = sigma[c(trt_ref, sel.dose[i])],
                              trtnames = c(trt_ref, sel.dose[i]))

    # Summarize the data
    n_int2_summary_data <- rbind(dat_int, dat_xtr) %>%
      dplyr::group_by(.data$Treatment) %>%
      dplyr::summarize(n = dplyr::n(), .groups = 'drop') %>%
      dplyr::mutate(sim = i)
    n.int2 <- dplyr::bind_rows(n.interim , n_int2_summary_data)







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
              th.fut = th.fut,
              th.eff = th.eff,
              th.prom = th.prom,
              gamma = gamma,
              nsim = nsim,
              method = method,
              n.interim = n.interim,
              n.final = n.final,
              interim = list(benefit = int.ben),
              simresults = simresults
  )

  class(out) <- "bcts"
  return(out)
}
