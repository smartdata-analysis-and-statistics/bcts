library(dplyr)

#' Conditional Power for Difference between Two Means
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


#' Title
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
eval_power <- function(n_int = 60, #20
                             n_pln = 150, #45 + 20, #
                             n_max = 180, #60 + 20, # 80
                              block_size = 6, # block size for randomization
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
                              nsim.ppos = 0, # no. of simulations to assess PPOS
                              nsim_int = 100, # no.of simulations at interim
                              num_chains = 4,
                              n.iter = 5000, #
                              n.adapt = 500,
                              perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                              progress.bar = "text", #
                              quiet = TRUE) {
  require(rjags)

  n_treat <- length(mu_t)
  n_arms <- n_treat + 1



  if (any(c(n_int, n_pln, n_max) < 0)) {
    stop("The sample size should be strictly greater than 0!")
  }
  if (n_int > n_pln) {
    stop("The sample size at interim cannot be greater than the planned sample size!")
  }

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


  # Keep track of observed parameter estimates at interim
  results_interim <- tibble(
    sim = integer(),
    Treatment = character(),
    trt = integer(),
    n_c = integer(),
    n_t = integer(),
    est_lower = double(),
    est = double(),
    est_upper = double(),
    est_se = double(),
    z_test = double(),
    I_k = double(),
    success = logical(),
    ppos.est = numeric(),
    N_fin = numeric(),
    fut.trt = logical()
  )

  # Assign the dynamic column names
  names(results_interim)[names(results_interim) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(results_interim)[names(results_interim) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))

  results_per_trial <- data.frame(sim = seq(nsim),
                                  n_int = NA,
                                  n_pln = NA,
                                  n_max = NA,
                                  fut.trig = FALSE,
                                  inc.ss = FALSE)


  # Keep track of the final sample size per arm
  n_final <- rep(n_pln, nsim)

  if (progress.bar == "text") {
    pb <- txtProgressBar(min = 0, max = nsim, initial = 0)
  }

  ##############################################################################
  # Run simulation
  # Instead of computing power for a fixed value of the parameters of interest,
  # in the Bayesian power computation, the parameters are sampled
  # (and thus change at each computation) from their assumed distribution.
  # For each sampled value, a unique simulated study is conducted, and the
  # success or failure of the study is assessed.
  ##############################################################################
  for (i in 1:nsim) {
    # Simulate all trial data
    dat_int <- generate_trial_data(N = n_int,
                               block_size = block_size,
                               mean = c(mu_c, mu_t),
                               sd = c(sd_c, sd_t),
                               trtnames = trt_names)

    results_per_trial$n_int[i] <- nrow(dat_int)
    results_per_trial$n_pln[i] <- n_pln
    results_per_trial$n_max[i] <- n_max

    # Assess success at the interim stage
    result_int <- eval_superiority(data = dat_int,
                                   margin = 0,
                                   gamma = gamma,
                                   num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)

    # Evaluate PPOS
    dat_post_int <- generate_trial_data(N = n_post_interim,
                                        block_size = block_size,
                                        mean = c(mu_c, mu_t),
                                        sd = c(sd_c, sd_t),
                                        trtnames = trt_names) %>% mutate(Y = NA)

    result_fin <- eval_superiority(rbind(dat_int, dat_post_int), margin = 0, gamma = gamma, num_chains = num_chains,
                     n.adapt = n.adapt, n.iter = n.iter, perc_burnin = perc_burnin)

    # Evaluate PPOS freq
    # Calculate pooled SD
    sigma <- sqrt(((20 - 1) * 1.10689 + (20 - 1) * 1.106016)/(40 - 2))

    ppos.seq1 <- 1/(2*sigma)
    ppos.seq2 <- sqrt(40/(130-40))
    ppos.seq3 <- (0.55950954*sqrt(130))-(2*sigma*qnorm(gamma))
    pnorm(ppos.seq1*ppos.seq2*ppos.seq3)



    # Store the interim results
    out_int <- cbind(sim = i, result_int$est %>% merge(ppos %>% select(Treatment, N_fin, ppos.est), by = "Treatment"))

    # Futility of individual treatments
    out_int <- out_int %>% mutate(fut.trt = ifelse(ppos.est < th.fut, TRUE, FALSE))
    results_per_trial$fut.trig[i] <- all(out_int$fut.trt) # Futility for all treatments?


    results_interim <- results_interim %>% add_row(out_int)

    # by default, the final sample size is the planned sample size





    if (progress.bar == "text") {
      setTxtProgressBar(pb, i)
    }
  }

  if (progress.bar == "text") {
    close(pb)
  }






  # Futility triggering
  mean(fut.trig)

  # Inc ss
  mean(inc.ss)

  # Final sample size
  mean(2*n_final + n_int)



  # Power at interim
  colMeans(sig_interim)






  # Calculate the probability that a certain dose is selected at interim
  P.sel <- table(sel_dose)/nsim
  names(P.sel) <- trt_names[-1]

  # Calculate average final sample size
  N.final <- colMeans(n_final)

  pop_mean <- c(mu_c, mu_t)
  pop_stdev <- c(sd_c, sd_t)
  names(pop_mean) <- names(pop_stdev) <- trt_names

  assumptions <- list(pop_mean = pop_mean,
                      pop_stdev = pop_stdev)
  result.start <- list(power = power_trial_start)
  result.interim <- list(P.fut.trig = mean(fut.trig),
                         P.sel = P.sel,
                         P.inc.ss = mean(inc.ss),
                         param_fut.trig = colMeans(param_mean_int[fut.trig,]),  # effect size at interim when futility is triggered
                         param_inc.ss = colMeans(param_mean_int[inc.ss,]), # effect size at interim when sample size increase is needed
                         psrf.est.0975 = apply(psrf_est_int, 2, function(x) quantile(x, 0.975)), #97.5 quantile of PSRF estimate
                         psrf.uci.0975 = apply(psrf_uci_int, 2, function(x) quantile(x, 0.975))) #97.5 quantile of PSRF upper CI

  result.final <- list(power = mean(hypothesis_testing),
                       N.final = N.final)

  out <- list(assumptions = assumptions,
              gamma = gamma,
              result.start = result.start,
              result.interim = result.interim,
              result.final = result.final)

  return(out)

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
#'
#' @return
#' @export
#'
#' @examples
eval_superiority <- function(data,
                             margin = 0,
                             gamma = 0.975, num_chains, n.adapt, n.iter, perc_burnin) {

  # Derive summary statistics at interim
  ad_int <- data  %>% group_by(Treatment) %>%
    summarize(trt = as.numeric(first(trt)),
              n = n(),
              mean = mean(Y),
              var = var(Y),
              tau_k = n/var(Y))

  # Get the maximum value of 'n'
  max_n <- max(ad_int$n)

  # Get the number of rows in 'ad_int'
  n_arms <- nrow(ad_int)

  # Create an empty dataframe 'Y' with the specified dimensions
  Y <- as.data.frame(matrix(NA, nrow = max_n, ncol = n_arms))

  # Assign column names based on the rows of 'ad_int' for clarity
  colnames(Y) <- ad_int$Treatment

  # Fill 'Y' with values from 'dat_int'
  for (treatment in levels(data$Treatment)) {
    treatment_data <- data %>% filter(Treatment == treatment)
    Y[1:nrow(treatment_data), treatment] <- treatment_data$Y
  }

  jags_data <- list(Y = Y,
                    n = ad_int %>% pull(n),
                    n_arms = n_arms,
                    z_alpha = qnorm(gamma))

  # Analyse the interim data
  jags_model_int <- jags.model(file = file.path("./inst/jags/eval_superiority_diff_means.jags"),
                               data = jags_data,
                               n.chains = num_chains,
                               n.adapt = n.adapt,
                               quiet = TRUE)

  ##Burnin stage
  update(jags_model_int, n.iter = ceiling(n.iter*perc_burnin), progress.bar = "none")

  ##Sampling after burnin
  fit_jags_int <- coda.samples(jags_model_int, variable.names =  c("trteff", "obs_diff", "R", "R_pred", "mu", "sigma"),
                               n.iter = floor(n.iter*(1 - perc_burnin)),
                               progress.bar = "none")

  ## Derive posterior distributions at interim
  psample_int <- as.data.frame(do.call(rbind, fit_jags_int))



  trt_est <- data.frame("Treatment" = character(),
                        "trt" = numeric(),
                        "n" = numeric(),
                        "est_lower" = numeric(),
                        "est" = numeric(),
                        "est_upper"  = numeric(),
                        "est_se" = numeric(), # Standard error of the mean
                        "z_test" = numeric(),
                        "I_k" = numeric(), # Information level
                        "pos" = numeric(),# Success based on parameter means (i.e., POS)
                        "ppos" = numeric())  # Trial success based on observed data (i.e.,)

  for (treat in 2:n_arms) {
    mudiff <- psample_int %>% pull(paste0("trteff[", treat, "]"))
    R <- psample_int %>% pull(paste0("R[", treat, "]"))
    R_PPOS <- psample_int %>% pull(paste0("R_pred[", treat, "]"))
    trtc <- (ad_int %>% filter(trt == treat) %>% pull(Treatment))
    trtn <- (ad_int %>% filter(trt == treat) %>% pull(trt))
    n_c <- ad_int %>% filter(trt == 1) %>% pull(n)
    n_t <- ad_int %>% filter(trt == trtn) %>% pull(n)


    trt_est <- trt_est %>% add_row(data.frame("Treatment" = trtc,
                                              "trt" = trtn,
                                              "n" = n_c + n_t,
                                              "est_lower" = quantile(mudiff, 0.025),
                                              "est" = mean(mudiff),
                                              "est_upper" = quantile(mudiff, 0.975),
                                              "est_se" = sd(mudiff),
                                              "z_test" = mean(mudiff)/sd(mudiff),
                                              "I_k" = 1/var(mudiff), # Should be similar to 1/((sigmasq_c/n_c)+(sigmasq_t/n_t))
                                              "pos" = mean(R), #quantile(mudiff, (1 - gamma)) > 0,
                                              "ppos" = mean(R_PPOS)))
  }

  # Assign the dynamic column names
  names(trt_est)[names(trt_est) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(trt_est)[names(trt_est) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))


  return(trt_est)
}



# Internal function to evaluate the PPoS at an interim stage
eval_ppos_mcmc <- function(data_interim,
                      n_post_interim, # For how many patients do we need to generate additional data?
                      gamma, num_chains, n.adapt, n.iter, perc_burnin,
                      nsim = 1000) {

  # Keep track of observed parameter estimates at interim
  results_final <- data.frame(sim = numeric(),
                              Treatment = character(),
                              trt = numeric(),
                              n_c = numeric(),
                              n_t = numeric(),
                              est_025 = numeric(),
                              est = numeric(),
                              est_975 = numeric(),
                              est_se = numeric(),
                              z_test = numeric(),
                              I_k = numeric(),
                              success = logical())


  for (j in seq(nsim)) {
    dat_post_int <- generate_trial_data(N = n_post_interim,
                                        block_size = block_size,
                                        mean = c(mu_c, mu_t),
                                        sd = c(sd_c, sd_t),
                                        trtnames = trt_names)
    dat_fin <- rbind(data_interim, dat_post_int)

    # Assess success at the final stage
    result_fin <- eval_superiority(data = dat_fin,
                                   margin = 0,
                                   gamma = gamma,
                                   num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)
    results_final <- results_final %>% add_row(cbind(sim = j, result_fin$est))
  }

  ppos <- results_final %>% group_by(Treatment) %>%
    summarize(n_fin = mean(n_c + n_t),
              ppos.est = mean(success),
              ppos.se = sqrt(ppos.est*(1 - ppos.est)/n()))

  return(ppos)
}

eval_ppos_eq <- function(est_int, # Results from the interim analysis
                         N_fin, # Number of patients in the trial at the final stage
                              gamma) {

  ppos <- data.frame(Treatment = character(),
                     n_int = numeric(),
                     N_fin = numeric(),
                     ppos.est = numeric())


  # Use approximation from Kundu et al. ignoring prior evidence

  n_c <- ds_int %>% filter(trt == 1) %>% summarize(n = n()) %>% pull(n)
  var_c <- var(ds_int %>% filter(trt == 1)%>% pull(Y))
  n_t <- ds_int %>% filter(trt != 1)%>% summarize(n = n()) %>% pull(n)
  var_t <- var(ds_int %>% filter(trt != 1)%>% pull(Y))

  delta <- mean(ds_int %>% filter(trt != 1)%>% pull(Y)) - mean(ds_int %>% filter(trt == 1)%>% pull(Y))

  a <- n_t/n_c # Allocation ratio for treated vs control
  r <- sqrt(((a + 1)**2)/a)

  n <- n_c + n_t
  N_c_post <- N_post/(a+1)
  N_t_post <- a*N_post/(a+1)
  N <- N_c_post + N_t_post + n_c + n_t

  # Calculate pooled SD
  sigma <- sqrt(((n_c - 1) * var_c + (n_t - 1) * var_t)/(n_c + n_t - 2))

  ppos.seq1 <- 1/(r*sigma)
  ppos.seq2 <- sqrt(n/(N-n))
  ppos.seq3 <- (delta*sqrt(N))-(r*sigma*qnorm(gamma))

  trtc <- (ds_int %>% filter(trt != 1)  %>% pull(Treatment))[1]

  data.frame(Treatment = trtc,
             n_int = n,
             N_fin = N,
             ppos.est = pnorm(ppos.seq1*ppos.seq2*ppos.seq3))


}


#' Title
#'
#' @param N Total number of patients needed
#' @param N_interim Total number of patients at interim
#' @param block_size Randomization block size
#' @param mean Effect size for each treatment
#' @param sd Standard deviation for each treatment
#'
#' @return A dataframe with simulated trial data
#' @export
#'
#' @examples
generate_trial_data <- function(N,
                                block_size, mean, sd, trtnames) {

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


  return(trial_data)
}

# Function to create a single block of random assignments
create_block <- function(block_size, treatment_groups) {
  assignments <- rep(treatment_groups, each = block_size / length(treatment_groups))
  sample(assignments)
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
