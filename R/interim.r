#' Apply Dose Selection and Sample Size Re-estimation at Interim
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
dose_selection_and_ssre <- function(dat_int, n_pln, n_max, trt_ref, trt_active,
                                    trt_rank,
                                    prioritize_low_rank = TRUE, gamma = 0.975,
                                    th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                                    method = "mcmc", num_chains = 4,
                                    n.iter = 5000, n.adapt = 500, perc_burnin = 0.2) {

  # Evaluate PPOS for each treatment at the planned sample size
  ppos <- setNames(rep(NA, length(trt_active)), trt_active)

  # Evaluate the effect size for each treatment
  benefit <- setNames(rep(NA, length(trt_active)), trt_active)

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
    benefit[trt_act_i] <- result_pln %>% pull("est")
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
                                       method = method, trt_ref = trt_ref)

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

  return(list(benefit = benefit, PPoS = ppos, result = result_final))
}


#' Apply Dose Selection at Interim
#'
#' @param dat_int Observed data at the interim stage
#' @param n_pln Planned sample size
#' @param trt_ref Character denoting the control treatment
#' @param trt_active Character vector denoting the active treatment names
#' @param trt_rank Optional named vector with preference ranking for each treatment (e.g., lower doses are preferred)
#' @param prioritize_low_rank  If multiple treatments have a posterior predictive power > 'th.eff', the treatment with lowest rank
#' will be selected at interim. Disabled by default
#' @param gamma Level to declare success. For a fixed design, gamma is typically chosen as 0.975 for a one-sided type-I error rate of 2.5%. However, an increase is usually needed because of dose selection.
#' @param th.fut Futility threshold for the predictive power calculated conditioned on the interim data
#' @param th.eff Efficacy threshold for the predictive power calculated conditioned on the interim data
#' @param method Method to estimate the power. Choose "bayes" for a fully Bayesian approach or "mcmc" for a frequentist approximation.
#' @param num_chains Number of MCMC chains
#' @param n.iter Number of MCMC iterations for estimation
#' @param n.adapt Number of MCMC iterations for adaptation
#' @param perc_burnin How many n.iter should be reserved for burnin?
#'
#' @return A list with PPoS and a dataframe containing information on the selected dose and final sample size
#'
#' @examples
#' # Generate interim data
#' mu <- c("Placebo" = 0, "Drug A" = 0.4, "Drug B" = 0.5)
#' sigma <- c("Placebo" = 1, "Drug A" = 1, "Drug B" = 1)
#' trt_ref <- "Placebo"
#' trt_active <- setdiff(names(mu), trt_ref)
#' dat_int <- sim_rct_normal(n = 20*3, mean = mu, sd = sigma)
#'
#' # Frequentist approach for dose selection
#' dose_selection(dat_int = dat_int, n_pln = 150, trt_ref = trt_ref,
#'                trt_active = trt_active, gamma = 0.98)
#'
#' # Bayesian approach for dose selection
#' dose_selection(dat_int = dat_int, n_pln = 150, trt_ref = trt_ref,
#'                trt_active = trt_active, method = "bayes")
#' @export
#'
#' @importFrom rlang .data
dose_selection <- function(dat_int, n_pln, trt_ref, trt_active,
                           trt_rank,
                           prioritize_low_rank = FALSE, gamma = 0.975,
                           th.fut = 0.2, th.eff = 0.9,
                           method = "mcmc", num_chains = 4,
                           n.iter = 5000, n.adapt = 500, perc_burnin = 0.2) {

  # Evaluate PPOS for each treatment at the planned sample size
  ppos <- setNames(rep(NA, length(trt_active)), trt_active)

  # Evaluate the effect size for each treatment
  benefit <- setNames(rep(NA, length(trt_active)), trt_active)

  # Evaluate the null hypothesis for each treatment
  rejectH0 <- setNames(rep(NA, length(trt_active)), trt_active)

  result_interim <- data.frame(fut.trig = NA,
                               inc.ss = FALSE,
                               sel.dose = "",
                               sel.dose.ppos = NA)

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
    benefit[trt_act_i] <- result_pln %>% pull("est")
    rejectH0[trt_act_i] <- result_pln %>% pull("rejectH0")
  }

  # Assess futility
  result_interim$fut.trig <- ifelse(max(ppos) < th.fut, TRUE, FALSE)

  if (max(ppos) < th.eff | prioritize_low_rank == FALSE | length(trt_active) == 1) {
    # Select the most effective dose
    result_interim$sel.dose <- names(ppos)[which.max(ppos)]
    result_interim$sel.dose.ppos <- max(ppos)
  } else {
    # Filter treatments that have ppos values above th.eff
    filtered_treatments <- ppos[ppos >= th.eff]

    # Get the treatment name with the lowest rank among the filtered treatments
    treatment_with_lowest_rank <- names(which.min(trt_rank[names(filtered_treatments)]))
    result_interim$sel.dose <- treatment_with_lowest_rank
    result_interim$sel.dose.ppos <- ppos[treatment_with_lowest_rank]
  }

  return(list(benefit = benefit, PPoS = ppos, rejectH0 = rejectH0, result = result_interim))
}


#' Apply Sample Size Re-estimation at Interim
#'
#' @param dat_int Observed data at the interim stage
#' @param n_pln Planned sample size
#' @param n_max Maximum sample size
#' @param trt_ref Character denoting the control treatment
#' @param trt_active Character denoting the active treatment
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
ssre <- function(dat_int, n_pln, n_max, trt_ref, trt_active,
                 trt_rank,
                 prioritize_low_rank = FALSE,
                 gamma = 0.975, th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                 method = "mcmc", num_chains = 4,
                 n.iter = 5000, n.adapt = 500, perc_burnin = 0.2) {

  result_interim <- data.frame(fut.trig = NA, inc.ss = NA, n.final = NA,
                               ppos_pln = NA, ppos_fin = NA, benefit = NA)

  # Simulate remaining trial data
  dat_xtr <- sim_rct_normal(n = n_pln - nrow(dat_int),
                            mean = rep(NA,2),
                            sd = rep(NA,2),
                            trtnames = c(trt_ref, trt_active))
  dat_pln <- rbind(dat_int %>% dplyr::filter(.data$Treatment %in% c(trt_ref, trt_active)), dat_xtr)
  dat_pln$Treatment <- droplevels(dat_pln$Treatment)

  result_pln <- eval_superiority(dat_pln, margin = 0, gamma = gamma,
                                 method = method, trt_ref = trt_ref,
                                 num_chains = num_chains,
                                 n.adapt = n.adapt, n.iter = n.iter,
                                 perc_burnin = perc_burnin)

  result_interim$ppos_pln <-  result_pln %>% pull("ppos")
  result_interim$benefit <- result_pln %>% pull("est")


  if (result_interim$ppos_pln < th.fut) {
    result_interim$fut.trig <- TRUE
    result_interim$inc.ss <- FALSE
    result_interim$ppos_fin <- result_interim$ppos_pln
    result_interim$n.final <- n_pln
  } else if (result_interim$ppos_pln >= th.prom & result_interim$ppos_pln < th.eff) {
    result_interim$fut.trig <- FALSE
    result_interim$inc.ss <- TRUE

    for (n_pln_new in (n_pln + 1):n_max) {
      dat_xtr <- sim_rct_normal(n = n_pln_new - nrow(dat_int),
                                mean = rep(NA, 2),
                                sd = rep(NA, 2),
                                trtnames = c(trt_ref, trt_active))

        dat_pln <- rbind(dat_int %>% dplyr::filter(.data$Treatment %in% c(trt_ref, trt_active)), dat_xtr)
        dat_pln$Treatment <- droplevels(dat_pln$Treatment)

        result_pln <- eval_superiority(dat_pln, margin = 0, gamma = gamma,
                                       method = method, trt_ref = trt_ref)

        ppos_eval <- result_pln %>% pull("ppos")

        result_interim$n.final <- n_pln_new
        result_interim$ppos_fin <- ppos_eval

        if (ppos_eval >= th.eff) {
          # Required sample size increase reached
          break;
        }
      }
    } else {
      result_interim$fut.trig <- FALSE
      result_interim$inc.ss <- FALSE
      result_interim$ppos_fin <- result_interim$ppos_pln
      result_interim$n.final <- n_pln
  }

  return(list(result = result_interim))
}
