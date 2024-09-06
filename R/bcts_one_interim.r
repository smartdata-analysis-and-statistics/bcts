#' Bayesian clinical trial simulation with One Interim Analysis (BCTS)
#'
#' @param n_int Interim sample size
#' @param n_pln Planned sample size
#' @param n_max Maximum sample size
#' @param mu Named vector with expected effect sizes
#' @param sigma Named vector with expected standard deviations
#' @param trt_ref Character denoting the control treatment
#' @param trt_rank Named vector with preference ranking for each treatment (e.g., lower doses are preferred) when multiple treatments have a posterior predictive power > 'th.eff'
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
#' @export
#' @return An object of class "bcts"
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr slice_head mutate select filter add_row pull
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom binom binom.confint
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
bcts_one_interim <- function(n_int, n_pln, n_max, mu, sigma, trt_ref, trt_rank = NULL,
                             gamma = 0.975,
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
                           n.final = NA,
                           est.final = NA,
                           est_lower.final = NA,
                           est_se.final = NA,
                           est_upper.final = NA,
                           rejectH0.final = NA)

  n.interim <- n.final <- data.frame(sim = numeric(), Treatment = character(), n = integer(), stringsAsFactors = FALSE)

  # Keep track of futility triggering (each column refers to an interim analysis)
  fut.trig <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = 1)), c("Int1"))

  # Keep track of benefit at interim
  int.ben <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = length(trt_active))), trt_active)

  # Keep track of PPoS for selected dose at interim
  int.PPoS <- setNames(data.frame(matrix(NA, nrow = nsim, ncol = 2)), c("Interim", "Final")) #interim 2

  # Keep track of sample size increase
  int.ss_inc <- rep(0, nsim)

  # Keep track of selected dose at first interim analysis
  sel.dose <- rep("", nsim)

  if (progress.bar == "text") {
    pb <- utils::txtProgressBar(min = 0, max = nsim, initial = 0)
  }

  for (i in seq(nsim)) {
    set.seed(seeds[i])

    # Simulate trial data at interim stage
    dat_int <- sim_rct_normal(n = n_int,
                              mean = mu,
                              sd = sigma,
                              trtnames = trt_names)

    ri <- dose_selection_and_ssre(dat_int = dat_int, n_pln = n_pln, n_max = n_max,
                                  trt_ref = trt_ref,
                                  trt_active = trt_active, trt_rank = trt_rank,
                                  gamma = gamma, th.fut = th.fut, th.eff = th.eff,
                                  th.prom = th.prom,
                                  method = method, num_chains = num_chains,
                                  n.iter = n.iter, n.adapt = n.adapt,
                                  perc_burnin = perc_burnin)

    fut.trig[i,1] <- ri$fut.trig # Futility triggering at interim 1
    sel.dose[i] <- ri$sel.dose # Selected dose at interim 1

    # Efficacy at interim
    int.ben[i,] <- ri$benefit

    # PPoS at interim
    int.PPoS[i,] <- ri$PPoS
    int.ss_inc[i] <- ri$N["Final"] - n_pln


    # Efficacy analysis in the final dataset
    dat_xtr <- sim_rct_normal(n = ri$N["Final"] - nrow(dat_int),
                              mean = mu[c(trt_ref, ri$sel.dose)],
                              sd = sigma[c(trt_ref, ri$sel.dose)],
                              trtnames = c(trt_ref, ri$sel.dose))
    dat_fin <- rbind(dat_int %>% dplyr::filter(.data$Treatment %in% c(trt_ref, ri$sel.dose)), dat_xtr)
    dat_fin$Treatment <- droplevels(dat_fin$Treatment)

    result_fin <- eval_superiority(data = dat_fin, margin = 0, gamma = gamma,
                                   trt_ref = trt_ref, method = method,
                                   num_chains = num_chains, n.adapt = n.adapt,
                                   n.iter = n.iter, perc_burnin = perc_burnin)

    simresults$seed[i] <- seeds[i]
    simresults$sel.dose[i] <- sel.dose[i]
    simresults$n.final[i] <- ri$N["Final"]
    simresults$fut.trig[i] <- any(fut.trig[i,])
    simresults$inc.ss[i] <- ri$inc.ss
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

  out <- list(no.looks = 1,
              trt_ref = trt_ref,
              trt_active = trt_active,
              mu = mu,
              sigma = sigma,
              th.fut = th.fut,
              th.eff = th.eff,
              th.prom = th.prom,
              gamma = gamma,
              nsim = nsim,
              method = method,
              fut.trig = fut.trig,
              sel.dose = sel.dose,
              PPos = list("Int1" = int.PPoS),
              benefit = list("Int1" = int.ben),
              ss_inc = list("Int1" = int.ss_inc),
              simresults = simresults
  )

  class(out) <- "bcts"
  return(out)
}
