library(dplyr)


source("R/simulation.r")

#' Estimate the Power for Difference between Two Means
#' The probability of success conditional on an assumed true treatment effect
#' @param n1
#' @param mu1
#' @param mu2
#' @param margin
#' @param sd1
#' @param sd2
#' @param kappa
#' @param alpha
#' @param alternative
#'
#' @description
#' For superiority, H0 is defined as mu1 - mu2 <= margin.
#'
#'
#' @return
#' @export
#'
#' @examples
pow_diff_means <- function(n1, # Sample size for the control arm
                             mu1, # Active treatment
                             mu2, # Control treatment
                             margin = 0, #non-inferiority/superiority margin
                             sd1,
                             sd2,
                             kappa = 1, # allocation ratio for the active treatment (n1/n2)
                             alpha = 0.05,
                             alternative = "superiority") {

  n2 <- ceiling(n1 / kappa)

  # calculate pooled SD
  sd_pooled <- sqrt(((n1 - 1)*sd1**2 + (n2 - 1)*sd2**2)/(n1 + n2 - 2))


  if (alternative == "superiority") {
    z_test <- (mu1 - mu2 - margin)/(sd_pooled * sqrt(1/n1 + 1/n2))
    z_alpha <- qnorm(1 - alpha)

    power <- pnorm(z_test - z_alpha)

    str_out <- "Difference between Two means\n"
    str_out <- paste0(str_out, "(Test for Noninferority/Superiority)\n")
    str_out <- paste0(str_out, "H0: mu1 - mu2 <= ", margin, "\n")
    str_out <- paste0(str_out, "H0: mu1 - mu2 > ", margin, "\n")
    str_out <- paste0(str_out, "------------------------------ \n")
    str_out <- paste0(str_out, " Statistical power = ", round(power,3), "\n")
    str_out <- paste0(str_out, " n1 = ", n1, "\n")
    str_out <- paste0(str_out, " n2 = ", n2, "\n")

    out <- list(power = power, txt = str_out)
    return(out)
  }
}


#' Title
#'
#' @param n_0
#' @param mu
#' @param sd
#' @param arm_names
#' @param delta
#' @param tar
#' @param alpha
#' @param power  Choose between "marginal", "disjunctive" or "conjunctive"
#' @param adjust alpha adjustment ( "bon","sid","no")
#' @param nsim
#'
#' @return
#' @export
#'
#' @examples
mcmc_power_diff_means <- function(n_0, # Sample size for the control arm
                             mu = c(0, 0, 0), # expected effect sizes
                             sd = c(1, 1, 1), # expected SD
                             arm_names = NULL,
                             delta = 0, #non-inferiority/superiority margin
                             tar = NULL, # Vector with allocation ratios
                             alpha = 0.05,
                             power = "marginal",
                             adjust = "bon", # multiple comparison correction (MCC)
                             nsim = 10000) {

  n_0 <- ifelse(n_0 < 1, 1, ceiling(n_0))
  K <- length(mu) - 1 # Number of active arms

  if (is.null(adjust)) {
    adjust <- "no"
  }
  if (is.null(arm_names)) {
    # Specify treatment names
    arm_names <- c("control", paste0("treatment_", seq(K)))
  }
  if (is.null(tar)) {
    tar <- rep(1, (K + 1))
  }
  n <- n_0 * tar
  names(n) <- arm_names

  # Determine the significance level
  if (K == 1) {
    if (adjust != "no") message("No correction for multiplicity is needed and thus will not be applied!")
    gamma <- alpha
  } else if (adjust == "bon") {
    gamma <- alpha/K
  } else if (adjust == "sid") {
    gamma <- 1 - (1 - alpha)^(1/K)
  } else {
    if (adjust == "none") warning("No correction for multiplicity was specified!")
    gamma <- alpha
  }

  p_k <- z_k <- sign <- matrix(data = NA, nrow = nsim, ncol = K)

  for (i in seq(nsim)) {
    y <- list()
    ybar <- rep(NA, (K + 1))
    for (k in seq(K + 1)) {
      y[[k]] <- rnorm(n = n[k], mean = mu[k], sd = sd[k])
      ybar[k] <- mean(y[[k]]) # Observed mean response

      # For each active treatment, make a comparison to the control treatment
      if (k > 1) {

        if (sd[k] == sd[1]) {
          I_k <- 1/((var(y[[1]])/n[1]) + (var(y[[k]])/n[k]))
          z_k[i, (k - 1)] <- (ybar[k] - ybar[1] - delta)*sqrt(I_k)
        } else {
          # Error variance of the observed response probabilities
          var_pooled <- (1/(n[1] + n[k] - 2)) * (sum((y[[1]] - ybar[1])**2) + sum((y[[k]] - ybar[k])**2))

          # Derive the z-score
          z_k[i, (k - 1)] <- (ybar[k] - ybar[1] - delta)/(sqrt(var_pooled)*sqrt(1/n[1] + 1/n[k]))
        }

        # Derive the p-value
        p_k[i, (k - 1)] <-  pnorm(z_k[i, (k - 1)], mean = 0, sd = 1, lower.tail = FALSE)

        # Determine whether to reject H0
        sign[i, (k - 1)] <- ifelse(z_k[i, (k - 1)] > qnorm(1 - gamma),1, 0)
      }
    }
  }

  # Check for each simulation how many hypotheses were rejected
  n_sign <-  rowSums(sign)
  n_sign_per_comparison <- colSums(sign)/nsim
  names(n_sign_per_comparison) <- arm_names[-1]

  # Estimate each power
  powers <- data.frame(matrix(NA, nrow = 1, ncol = K + 2))
  colnames(powers) <- c("disjunctive", "conjunctive", paste0("marginal_", arm_names[-1]))

  powers$disjunctive <- sum(n_sign > 0)/nsim
  powers$conjunctive = sum(n_sign == K)/nsim

  act_trt <- arm_names[-1]

  for (k in 1:K) {
    powers[1,paste0("marginal_", act_trt[k])] <- n_sign_per_comparison[k]
  }

  # Case evaluation to assign the value of power
  est.power <- case_when(
    power == "marginal" ~ min(n_sign_per_comparison),
    power == "disjunctive" ~ powers$disjunctive,
    power == "conjunctive" ~ powers$conjunctive,
    TRUE ~ NA_real_ # Default case if none of the above matches
  )

  dat <- data.frame(z_k)
  names(dat) <- paste0("z_",arm_names[-1])

  out <- list(power = est.power,
              powers = powers,
              n_0 = n_0,
              n = n,
              N = sum(n),
              K = K, # Number of active treatments
              alpha = alpha,
              power_type = power,
              p_crit = gamma,
              adjust = adjust,
              z_test_lower = qnorm(1 - gamma))
  return(out)
}



n_diff_means <- function(mu, sd, delta, alpha, beta, power_type, adjust) {
  uniroot(function(x) (mcmc_power_diff_means(n_0 = x,
                                        mu = mu,
                                        sd = sd,
                                        delta = delta,
                                        adjust = adjust,
                                        power = power_type,
                                        alpha = alpha)),
          beta = beta)
}


#' Bayesian clinical trial simulation (BCTS)
#'
#' @param n_int Planned sample size at interim (control arm)
#' @param n_0_pln Planned sample size for the control group
#' @param n_0_max Maximum sample size for the control group after sample size increase
#' @param mu_c Expected effect size for the control treatment
#' @param sd_c Expected standard deviation for the control treatment
#' @param mu_t Named vector with expected effect sizes for the active treatments
#' @param sd_t Named vector with expected standard deviations for the active treatments
#' @param prioritize_low_trteff  If multiple treatments have a posterior predictive power > 'th.eff', the treatment with lowest effect size
#' will be selected at interim  (prioritize_low_trteff = TRUE).
#' @param gamma Level to declare success. For a fixed design, gamma is typically chosen as 0.975 for a one-sided type-I error rate of 2.5%. However, an increase is usually needed because of dose selection.
#' @param th.fut Futility threshold for the predictive power calculated conditioned on the interim data
#' @param th.eff Efficacy threshold for the predictive power calculated conditioned on the interim data
#' @param th.prom Predictive power threshold that would trigger sample size increase
#' @param nsim Number of trials to simulate (default: 1000)
#' @param nsim.ppos  Number of simulations to assess PPOS (set to 0 to use approximation)
#' @param num_chains Number of MCMC chains
#' @param n.iter Number of MCMC iterations for estimation
#' @param n.adapt Number of MCMC iterations for adaptation
#' @param perc_burnin How many n.iter should be reserved for burnin?
#' @param progress.bar Type of progress bar. Possible values are "text", "gui", and "none"
#' @param quiet Display messaging output
#'
#' @description
#' The null hypothesis can be rejected when the posterior probability that $\mu_t$ - $\mu_c$ exceeds a high probability threshold,
#' i.e. $Pr(\mu_t - \mu_c > 0|D) > $\gamm$)$.
#'
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' mu_t = c("Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5)
#' sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1)
#'
#' eval_power(n_int = 60, n_pln = 150, n_max = 180, block_size = 6,
#'            mu_c = 0, sd_c = 1, mu_t = mu_t, sd_t = sd_t, gamma = 0.98, nsim.ppos = 0)
#' }
bcts <- function(n_per_arm_int = 20,
                             n_per_arm_pln = 65,
                             n_per_arm_max = 80,
                       block.sizes = 1, # block size for randomization
                              mu_c = 0,
                              sd_c = 1,
                              mu_t = c("Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
                              sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
                              prioritize_low_trteff = TRUE,
                              gamma = 0.98, #   %
                              th.fut = 0.2, ##
                              th.eff = 0.9, ##
                              th.prom = 0.5, ##
                              nsim = 1000, # no. of simulatoins at trial start
                              num_chains = 4,
                              n.iter = 5000, #
                              n.adapt = 500,
                              perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                              progress.bar = "text", #
                              quiet = TRUE) {
  require(rjags)

  # Sample the random seeds
  seeds <- sample(seq(nsim*10,nsim*20), size = nsim)


  n_treat <- length(mu_t)
  n_arms <- n_treat + 1

  # Calculate the sample sizes
  n_pln <- 2*n_per_arm_pln + (length(mu_t) - 1)*n_per_arm_int
  n_max <- 2*n_per_arm_max + (length(mu_t) - 1)*n_per_arm_int


  # Check if all elements are equal
  if (all(mu_t == mu_c) & !quiet) {
    message("Evaluating the type-I error")
  }

  # Extract the treatment names
  if (is.null(names(mu_t))) {
    names(mu_t) <- names(sd_t) <- paste0("trt", seq(mu_t))
  }
  trt_names <- c("Control", names(mu_t))
  mu <- c(mu_c, mu_t)
  sigma <- c(sd_c, sd_t)
  names(mu) <- names(sigma) <- trt_names

  simresults <- data.frame(sim = seq(nsim),
                           seed = NA,
                           fut.trig = FALSE,
                           inc.ss = FALSE,
                           sel.dose = NA,
                           sel.dose.ppos = NA,
                           n_per_arm_fin = NA,
                           sig_fin = NA)

  results_interim <- results_final <- results_max <- data.frame(sim = numeric(),
                           Treatment = character(),
                           trt = numeric(),
                           n = numeric(),
                           N = numeric(),
                           est_lower = numeric(),
                           est = numeric(),
                           est_upper = numeric(),
                           est_se = numeric(),
                           zb = numeric(),
                           zf = numeric(),
                           Ib = numeric(),
                           If = numeric(),
                           sigb = numeric(),
                           sigf = numeric(),
                           pposb_mcmc = numeric(),
                           pposb_approx = numeric())
  names(results_interim)[names(results_interim) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(results_interim)[names(results_interim) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))
  names(results_final)[names(results_final) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(results_final)[names(results_final) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))
  names(results_max)[names(results_max) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(results_max)[names(results_max) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))


  if (progress.bar == "text") {
    pb <- txtProgressBar(min = 0, max = nsim, initial = 0)
  }

  for (i in 1:nsim) {
    simresults$seed[i] <- seeds[i]
    set.seed(seeds[i])

    # Simulate all trial data for the maximum sample size
    dat_max <- sim_rct_normal(n = n_per_arm_max*n_arms,
                              mean = c(mu_c, mu_t),
                              sd = c(sd_c, sd_t),
                              trtnames = trt_names,
                              block.sizes = block.sizes)

    result_max <- eval_superiority(dat_max,
                                   margin = 0, # Superiority margin
                                   gamma = gamma, num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)

    results_max <- results_max %>% add_row(cbind(sim = i, result_max))

    # Evaluate PPOS at the planned sample size
    dat_pln <- dat_max %>% slice_head(n = n_per_arm_pln*n_arms) %>%
      mutate(Y = ifelse(id > n_per_arm_int*n_arms, NA, Y))

    result_pln <- eval_superiority(dat_pln,
                                   margin = 0, # Superiority margin
                                   gamma = gamma, num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)

    results_interim <- results_interim %>% add_row(cbind(sim = i, result_pln))




    if (max(result_pln$pposb_mcmc) < th.fut) {
      simresults$fut.trig[i] <- TRUE
      simresults$inc.ss[i] <- FALSE
      simresults$n_per_arm_fin[i] <- n_per_arm_pln
      simresults$sel.dose[i] <- result_pln$trt[order(result_pln$pposb_mcmc, decreasing = TRUE)[1]]
      simresults$sel.dose.ppos[i] <- result_pln$pposb_mcmc[order(result_pln$pposb_mcmc, decreasing = TRUE)[1]]
    } else if (max(result_pln$pposb_mcmc) >= th.prom & max(result_pln$pposb_mcmc) < th.eff) {
      simresults$fut.trig[i] <- FALSE
      simresults$inc.ss[i] <- TRUE

      # Evaluate predictive power for each sample size
      ppos_by_n <- data.frame(n_per_arm_new = numeric(),
                              Treatment = character(),
                              trt = numeric(),
                              n = numeric(),
                              N = numeric(),
                              pposb_mcmc = numeric(),
                              pposb_approx = numeric())

      for (n_per_arm_new in (n_per_arm_pln + 1):n_per_arm_max) {

        dat_pln_new <- dat_max %>% slice_head(n = n_per_arm_new*n_arms) %>%
          mutate(Y = ifelse(id > n_per_arm_int*n_arms, NA, Y))

        result_pln_new <- eval_superiority(dat_pln_new,
                                           margin = 0, # Superiority margin
                                           gamma = gamma, num_chains = num_chains,
                                           n.adapt = n.adapt, n.iter = n.iter,
                                           perc_burnin = perc_burnin)

        ppos_by_n <- ppos_by_n %>% add_row(cbind(n_per_arm_new, result_pln_new %>% select(Treatment, trt, n, N, pposb_mcmc, pposb_approx)))
      }

      # Select the lowest sample size resulting in a PPOS >= th.eff
      filtered_ppos <- subset(ppos_by_n, pposb_mcmc >= th.eff)

      if (nrow(filtered_ppos) > 0) {
        sel_row <- filtered_ppos[which.min(filtered_ppos$n_per_arm_new),]
        simresults$sel.dose[i] <- sel_row$trt
        simresults$sel.dose.ppos[i] <- sel_row$pposb_mcmc
        simresults$n_per_arm_fin[i] <- sel_row$n_per_arm_new
      } else {
        # Increasing the sample size did not help to attain target power
        # Keep the maximum sample size
        simresults$sel.dose[i] <- result_pln$trt[order(result_pln$pposb_mcmc, decreasing = TRUE)[1]]
        simresults$sel.dose.ppos[i] <- result_pln$pposb_mcmc[order(result_pln$pposb_mcmc, decreasing = TRUE)[1]]
        simresults$n_per_arm_fin[i] <- n_per_arm_max
      }
    } else {
      simresults$fut.trig[i] <- FALSE
      simresults$inc.ss[i] <- FALSE
      simresults$n_per_arm_fin[i] <- n_per_arm_pln

      # If both predictive powers exceed 0.9 then
      # choose Low dose otherwise choose the dose with higher
      # predictive power
      if (sum((result_pln$pposb_mcmc >= th.eff)) >= 2 & prioritize_low_trteff) {
        simresults$sel.dose[i] <- result_pln$trt[order(mu_t[which(result_pln$pposb_mcmc >= th.eff)], decreasing = FALSE)[1]]
        simresults$sel.dose.ppos[i] <- result_pln$pposb_mcmc[order(mu_t[which(result_pln$pposb_mcmc >= th.eff)], decreasing = FALSE)[1]]
      } else {
        simresults$sel.dose[i] <- result_pln$trt[order(result_pln$pposb_mcmc, decreasing = TRUE)[1]]
        simresults$sel.dose.ppos[i] <- result_pln$pposb_mcmc[order(result_pln$pposb_mcmc, decreasing = TRUE)[1]]
      }
    }

    # Keep the reference treatment and the selected treatment
    dat_fin <- dat_max %>% filter(trt %in% c(1,simresults$sel.dose[i])) %>%
      slice_head(n = simresults$n_per_arm_fin[i]*2)

    result_fin <- eval_superiority(dat_fin,
                                   margin = 0, # Superiority margin
                                   gamma = gamma, num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)

    results_final <- results_final %>% add_row(cbind(sim = i, result_fin))

    simresults$sig_fin[i] <- result_fin %>% pull(sigb)

    if (progress.bar == "text") {
      setTxtProgressBar(pb, i)
    }
  }

  if (progress.bar == "text") {
    close(pb)
  }

  out.assumptions <- list(mu_c = mu_c, mu_t = mu_t, sd_c = sd_c, sd_t = sd_t)
  out.n <- list(interim = n_per_arm_int*n_arms)

  out <- list(assumptions = out.assumptions,
              n = out.n,
              gamma = gamma,
              power = mean(simresults$sig),
              simresults = simresults,
              results_interim = results_interim,
              results_max = results_max,
              results_final = results_final
              )

  class(out) <- "ssbayes"
  return(out)
}

eval_power <- function(dat, ...) {

}


#' Evaluate the posterior probability of treatment superiority
#'
#' @param data Trial data
#' @param margin Superiority margin (default: 0)
#' @param gamma Level to declare success (default: 0.975)
#' @param num_chains
#' @param n.adapt
#' @param n.iter
#' @param perc_burnin
#'
#' @description
#' Superiority is established if Pr(mean(treat)-mean(control)>margin) > gamma
#'
#' @references
#' O’Hagan A, Stevens JW, Campbell MJ. Assurance in clinical trial design. Pharmaceut Statist. 2005 Jul;4(3):187–201.
#'
#'
#' @return
#' @export
#'
#' @examples
eval_superiority <- function(data,
                             margin = 0,
                             gamma, #0.975
                             num_chains,
                             n.adapt, n.iter, perc_burnin) {

  treatments <- unique(data$trt)
  z_crit <- qnorm(gamma)



  jags.config <- prepare_jags_ppos(dat = data, gamma = gamma)


  # Analyse the interim data
  jags_model <- jags.model(file = textConnection(jags.config$model_string),
                           data = jags.config$jags_data,
                           n.chains = num_chains,
                           n.adapt = n.adapt,
                           quiet = TRUE)

  ##Burnin stage
  update(jags_model, n.iter = ceiling(n.iter*perc_burnin), progress.bar = "none")

  ##Sampling after burnin
  fit_jags_int <- coda.samples(jags_model, variable.names =  jags.config$variable.names,
                               n.iter = floor(n.iter*(1 - perc_burnin)),
                               progress.bar = "none")

  ## Derive posterior distributions at interim
  psample_int <- as.data.frame(do.call(rbind, fit_jags_int))
  sigma_c <- mean(psample_int %>% pull(paste0("sigma[1]")))

  trt_est <- data.frame("Treatment" = character(),
                        "trt" = numeric(),
                        "n" = numeric(),
                        "N" = numeric(),
                        "est_lower" = numeric(),
                        "est" = numeric(),
                        "est_upper"  = numeric(),
                        "est_se" = numeric(), # Standard error of the mean
                        "zb" = numeric(),
                        "zf" = numeric(),
                        "Ib" = numeric(), # Information level
                        "If" = numeric(), # Information level
                        "sigb" = numeric(), # Significance level of "est"
                        "sigf" = numeric(), # Significance level of "est"
                        "pposb_mcmc" = numeric(),
                        "pposb_approx" = numeric())  # Trial success based on observed data (i.e.,)

  for (treat in treatments[treatments != 1]) {
    mudiff <- psample_int %>% pull(paste0("trteff[", treat, "]"))
    sigma_t <- mean(psample_int %>% pull(paste0("sigma[", treat, "]")))
    PPOS <- psample_int %>% pull(paste0("ppos[", treat, "]"))
    trtc <-  (data %>% filter(trt == treat) %>% pull(Treatment))[1]
    n <-  nrow(data %>% filter(trt %in% c(1, treat) & !is.na(Y)))
    N <-  nrow(data %>% filter(trt %in% c(1, treat)))

    # Evaluate significance using frequentist method
    Yc <- data %>% filter(trt == 1 & !is.na(Y)) %>% pull(Y)
    Yt <- data %>% filter(trt == treat & !is.na(Y)) %>% pull(Y)
    sepooled <- sqrt(var(Yc)/length(Yc) + var(Yt)/length(Yt))
    sdpooled <- sqrt((sigma_t**2 + sigma_c**2)/2)
    z_freq <- (mean(Yt) - mean(Yc))/sepooled

    alr <- length(Yt)/length(Yc) # Allocation ratio
    r <- sqrt(((alr + 1)^2) / alr)

    # Calculate Predictive Power of success (PPoS)
    pposb_approx <- pnorm((1/(r*sdpooled))*sqrt(n/(N-n))* ((mean(mudiff)-margin)*sqrt(N)-r*sdpooled*z_crit))

    trt_est <- trt_est %>% add_row(data.frame("Treatment" = trtc,
                                              "trt" = as.numeric(treat),
                                              "n" = n,
                                              "N" = N,
                                              "est_lower" = quantile(mudiff, 1 - gamma),
                                              "est" = mean(mudiff),
                                              "est_upper" = quantile(mudiff, gamma),
                                              "est_se" = sd(mudiff),
                                              "zb" = mean(mudiff)/sd(mudiff),
                                              "zf" = z_freq,
                                              "Ib" = 1/var(mudiff),
                                              "If" = 1/(sepooled**2),
                                              "sigb" = quantile(mudiff, (1 - gamma)) > 0,
                                              "sigf" = z_freq > z_crit,
                                              "pposb_mcmc" = mean(PPOS),
                                              "pposb_approx" = pposb_approx))
  }

  # Assign the dynamic column names
  names(trt_est)[names(trt_est) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(trt_est)[names(trt_est) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))


  return(trt_est)
}

#' Assess the Predictive Power of success (PPoS)in an empirical manner
#'
#' @param interim_data
#' @param N
#' @param mu_c
#' @param mu_t
#' @param sd_c
#' @param sd_t
#' @param block.sizes
#' @param margin
#' @param gamma
#' @param num_chains
#' @param n.adapt
#' @param n.iter
#' @param perc_burnin
#'
#' @return
#' @export
#'
#' @examples
eval_predictive_power <- function(interim_data,
                                  N_per_arm,  # Planned sample size per arm
                                  mu_c, mu_t, sd_c, sd_t,
                                  block.sizes, margin = 0,
                                  gamma = 0.975,
                                  num_chains,
                                  n.adapt, n.iter, perc_burnin) {

  n_arms <- length(unique(interim_data$trt))

  sign <- matrix(NA, nrow = 1000, ncol = n_arms)


  for (i in 1 : 1000) {
    # Generate extra data
    ds <- sim_rct_normal(n = N_per_arm*n_arms - nrow(interim_data),
                         mean = c(mu_c, mu_t),
                         sd = c(sd_c, sd_t),
                         trtnames = trt_names,
                         block.sizes = block.sizes)
    rs <- eval_superiority(rbind(interim_data, ds),
                          margin = 0, # Superiority margin
                          gamma = gamma, num_chains = num_chains,
                          n.adapt = n.adapt, n.iter = n.iter,
                          perc_burnin = perc_burnin)

    for (k in 2:n_arms) {
      sign[i,rs$trt] = rs$sigb
    }
  }

  return(sign)
}




#' Title
#'
#' @param N Total number of patients needed
#' @param block_size Randomization block size
#' @param mean Effect size for each treatment
#' @param sd Standard deviation for each treatment
#' @param trtnames Treatment names
#' @param censor_Y Do we observe the outcomes?
#'
#' @return A dataframe with simulated trial data
#' @export
#'
#' @examples
generate_trial_data <- function(N,
                                block_size,
                                mean,
                                sd,
                                trtnames,
                                censor_Y = FALSE) {

  # Function to create a single block of random assignments
  create_block <- function(block_size, treatment_groups) {
    assignments <- rep(treatment_groups, each = block_size / length(treatment_groups))
    sample(assignments)
  }

  n_arms <- length(mean)

  # Number of randomization blocks needed until interim
  num_blocks_int <- ceiling(N/block_size)

  if (block_size %% n_arms != 0) {
    stop("Please adjust the block size to ensure it is a multiple of the number of arms")
  }

  # Generate the full randomization sequence
  randomization_sequence <- unlist(lapply(seq(num_blocks_int), function(x) create_block(block_size, trtnames)))

  # Create a data frame to store the simulated trial data
  trial_data <- data.frame(
    Participant = seq(N),
    Treatment = randomization_sequence[1:N],
    trt = NA,
    Y = NA
  )

  # Count number of patients per treatment
  n_per_arm <- trial_data %>% group_by(Treatment) %>% summarize(n = n())

  for (treat in seq(mean)) {
    Y_sim <- rnorm(n_per_arm %>% filter(Treatment == trtnames[treat]) %>% pull(n), mean[treat], sd[treat])
    trial_data <- trial_data %>% mutate(Y = ifelse(Treatment == trtnames[treat], Y_sim, Y))
  }
  trial_data <- trial_data %>% mutate(Treatment = factor(Treatment),
                                      trt = factor(Treatment,
                                                         labels = seq(trtnames),
                                                         levels = trtnames))

  if (censor_Y) {
    trial_data <- trial_data %>% mutate(Y = NA)
  }

  return(trial_data)
}





bpower_diff_means <- function(n_c,
                              n_t,
                              mu_c,
                              mu_t,
                              dgen_priorsd_mu_c = 0, # Data generation prior for mu_c (choose 0 for a frequentist approach)
                              dgen_priorsd_mu_t = 0,
                              sd_c,
                              sd_t,
                              gamma = 0.975, # Probability needed to declare success
                              nsim = 1000, #
                              num_chains = 4,
                              n.iter = 5000, #
                              n.adapt = 500,
                              perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                              progress.bar = "text", #
                              model_dir = ".") {

  require(rjags)


  pb <- txtProgressBar(min = 0, max = nsim, initial = 0)
  trteff_pos <- rep(NA, nsim)

  # Sample the "True" treatment effect for each simulation
  sim_mu_c <- rnorm(nsim, mean = mu_c, sd = dgen_priorsd_mu_c)
  sim_mu_t <- rnorm(nsim, mean = mu_t, sd = dgen_priorsd_mu_t)

  for (i in 1:nsim) {
    Y_c <- rnorm(n_c, mean = sim_mu_c[i], sd = sd_c)
    Y_t <- rnorm(n_t, mean = sim_mu_t[i], sd = sd_t)

    # Estimate error variance
    evar_c <- var(Y_c)/n_c
    evar_t <- var(Y_t)/n_t

    jags_data <- list(mean_c = mean(Y_c),
                      mean_t = mean(Y_t),
                      tau_c = 1/evar_c,
                      tau_t = 1/evar_t)

    # Parameters to monitor
    params <- c("trteff")

    # Analyse the interim data
    jags_model_int <- jags.model(file = file.path(model_dir, "inst/jags/diff_means_twoarm.jags"),
                                 data = jags_data,
                                 n.chains = num_chains,
                                 n.adapt = n.adapt,
                                 quiet = TRUE)

    ##Burnin stage
    update(jags_model_int, n.iter = ceiling(n.iter*perc_burnin), progress.bar = "none")

    ##Sampling after burnin
    fit_jags_int <- coda.samples(jags_model_int, params,
                                 n.iter = floor(n.iter*(1 - perc_burnin)),
                                 progress.bar = "none")

    ## Derive posterior distributions at interim
    psample_int <- as.data.frame(do.call(rbind, fit_jags_int))

    trteff_pos[i] <- quantile(psample_int %>% pull("trteff"), 1 - gamma) > 0
    setTxtProgressBar(pb,i)
  }
  close(pb)

  PoS <- mean(trteff_pos)

  return(PoS)
}

prepare_jags_ppos <- function(dat, gamma) {
  model_string <- "model {\n"

  treatments <-  unique(dat$trt)

  for (treat in treatments) {
    model_string <- paste0(model_string, "   for (i in 1:length(Y", treat, ")) {\n")
    model_string <- paste0(model_string, "      Y", treat, "[i] ~ dnorm(mu[", treat, "], prec[", treat, "])\n")
    model_string <- paste0(model_string, "   }\n")
    model_string <- paste0(model_string, "   mu[", treat,"] ~ dnorm(0, 0.001)\n")
    model_string <- paste0(model_string, "   prec[", treat,"] ~ dgamma(0.001, 0.001)\n")
    model_string <- paste0(model_string, "   sigma[", treat,"] <- 1 / sqrt(prec[", treat,"])\n\n")

    if (treat > 1) {
      model_string <- paste0(model_string, "   trteff[", treat, "] <- mu[", treat, "] - mu[1]\n")
      model_string <- paste0(model_string, "   tau[", treat, "] <- sqrt(sigma[", treat, "]**2/length(Y", treat, ") + sigma[1]**2/length(Y1))\n")
      model_string <- paste0(model_string, "   ppos[", treat, "] <- trteff[", treat, "] > tau[", treat, "]*z_alpha\n\n")
    }
  }

  # See https://www.math.chalmers.se/Stat/Grundutb/GU/MSA620/S20/Assurance.pdf

  model_string <- paste0(model_string, "}")

  # Split the data frame by 'trt' and extract the 'Y' values
  Y_list <- split(dat$Y, dat$trt)

  # Rename the list elements to Y1, Y2, Y3, etc.
  names(Y_list) <- paste0("Y", names(Y_list))

  # Remove empty elements
  Y_list <- Y_list[sapply(Y_list, length) > 0]

  Y_list$z_alpha <- qnorm(gamma)

  variable.names =  c("trteff", "mu", "sigma", "ppos")


  return(list(model_string = model_string, jags_data = Y_list,
              variable.names = variable.names))
  }



