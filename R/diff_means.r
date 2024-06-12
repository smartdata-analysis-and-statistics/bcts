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
cp_diff_means <- function(n1, # Sample size for the control arm
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

print.estsample <- function(x, ...) {

  cat(paste("Design parameters for a 1 stage trial with ",
            x$K, " active treatments\n\n", sep = ""))

  cat(paste("Total sample size:", x$N, "\n"))


  if (x$correction == "bonferroni") {
    descrip.gamma <- paste0("The maximum probability of incorrectly rejecting at least one of the ", x$K, " null hypotheses is at most ", x$alpha, ".")
  }
  cat(descrip.gamma)

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
bpower_adapt_diff_means <- function(n_int = 20, # 16
                                    n_pln = 75, # NOte that one arm should be dropped
                              n_max = 90, # 80
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
                              nsim_int = 100, # no.of simulations at interim
                              num_chains = 4,
                              n.iter = 5000, #
                              n.adapt = 500,
                              perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                              progress.bar = "text", #
                              quiet = TRUE) {
  require(rjags)

  ntreat <- length(mu_t)
  n_arms <- ntreat + 1



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

  power_trial_start <- calc_cp(mu_c = mu_c, mu_t = mu_t,
                                  sd_c = sd_c, sd_t = sd_t, margin = 0,
                                  n_int_c = n_int,
                                  n_pln_c = n_pln,
                                  n_max_c = n_max,
                                  kappa = 1, alpha = 0.025)




  ##############################################################################
  # Simulation results to monitor
  ##############################################################################

  # Keep track of success at interim
  sig_interim <- matrix(nrow = nsim, ncol = ntreat)
  colnames(sig_interim) <-  trt_names[-1]

  # Number of sample sizes to evaluate
  sig_pln <- lapply(1:length(n_pln:n_max), function(x) matrix(NA, nrow = nsim, ncol = 3))

  # Keep track of PPoS at planned sample size
  PPoS_pln <- matrix(nrow = nsim, ncol = ntreat)
  colnames(PPoS_pln) <- trt_names[-1]

  # keep track of predictive power for each dose at interim
  pred_power_int <- matrix(nrow = nsim, ncol = ntreat)
  colnames(pred_power_int) <- trt_names[-1]

  # keep track of Potential scale reduction factor at interim
  psrf_est_int <- psrf_uci_int <- matrix(nrow = nsim, ncol = (1 + 2*ntreat))
  colnames(psrf_est_int) <- colnames(psrf_uci_int) <- c("mu_c",
                              paste0("mu_t[",seq(ntreat), "]"),
                              paste0("trteff[",seq(ntreat), "]"))

  # Keep track of observed parameter estimates at interim
  param_mean_int <- matrix(nrow = nsim, ncol = (1 + 2*ntreat))
  colnames(param_mean_int) <- colnames(psrf_est_int)

  # keep track of futility triggering at interim
  fut.trig <- rep(FALSE, nsim)

  # Keep track of the selected dose at interim
  sel_dose <- rep(NA, nsim)

  # Keep track of whether sample size is increased
  inc.ss <- rep(FALSE, nsim)



  # Keep track of hypothesis test at final stage
  hypothesis_testing <- rep(NA, nsim)

  # Keep track of the final sample size of control and selected treatment
  n_final <- matrix(nrow = nsim, ncol = 2)
  colnames(n_final) <- c("control", "selected arm")

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
    dat_fin <- generate_trial_data(N = n_max * n_arms,
                               block_size = block_size,
                               mean = c(mu_c, mu_t),
                               sd = c(sd_c, sd_t),
                               trtnames = trt_names)

    # Set observations after interim to NA
    dat_fin <- dat_fin %>% mutate(Y = ifelse(Participant > n_int*n_arms, NA, Y))

    # Select the interim data
    dat_int <- dat_fin %>% filter(Participant <= n_int*n_arms)

    # Assess success at the interim stage
    pos_int <- pos_interim_analysis(data = dat_int, gamma = gamma,
                                    num_chains = num_chains,
                                    n.adapt = n.adapt, n.iter = n.iter,
                                    perc_burnin = perc_burnin)
    sig_interim[i,names(pos_int)] <- pos_int

    # Assess success at trial end for all possible sample sizes
    ppos_int <- predict_success(data = dat_fin,
                                N_planned = n_pln, # Number of planned subjects per arm
                                gamma = gamma,
                                num_chains = num_chains,
                                n.adapt = n.adapt, n.iter = n.iter,
                                perc_burnin = perc_burnin)
    for (ni in seq(n_pln:n_max)) {
      sig_pln[[ni]][i,1] <- ppos_int[ni,1] # Sample size
      sig_pln[[ni]][i,2] <- ppos_int[ni,2] # 15 mg
      sig_pln[[ni]][i,3] <- ppos_int[ni,3] # 30 mg
    }

    if (progress.bar == "text") {
      setTxtProgressBar(pb, i)
    }
  }

  if (progress.bar == "text") {
    close(pb)
  }

  # Power at interim
  colMeans(sig_interim)

  # Predictive power at interim
  for (ni in seq(n_pln:n_max)) {
    sig_pln[[ni]][i,1]
    mean(sig_pln[[ni]][i,2])
    mean(sig_pln[[ni]][i,3])
  }

  ## TODO: run trial at new sample size


  colMeans(sig_pln)






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

pos_interim_analysis <- function(data, gamma = 0.975, num_chains, n.adapt, n.iter, perc_burnin) {

  n_arms <- length(unique(data$trt))
  n_treat <- n_arms - 1

  # Derive summary statistics at interim
  ad_int <- data  %>% group_by(Treatment) %>%
    summarize(n = n(),
              mean = mean(Y), var = var(Y),
              tau = n/var(Y))

  jags_data <- list(mean_c = ad_int %>% filter(Treatment == "Control") %>% pull(mean),
                    tau_c = ad_int %>% filter(Treatment == "Control") %>% pull(tau),
                    mean_t = ad_int %>% filter(Treatment != "Control") %>% pull(mean),
                    tau_t = ad_int %>% filter(Treatment != "Control") %>% pull(tau))

  # Parameters to monitor
  params <- c("trteff")

  # Analyse the interim data
  jags_model_int <- jags.model(file = file.path("./inst/jags/diff_means_multiarm.jags"),
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

  ## Derive potential scale reduction factor
  gd_int <- gelman.diag(fit_jags_int, multivariate = F)
  psrf.est <-  gd_int$psrf[,"Point est."] # Point estimate
  psrf.uci <-  gd_int$psrf[,"Upper C.I."] # Upper C.I.

  ## Derive posterior distributions at interim
  psample_int <- as.data.frame(do.call(rbind, fit_jags_int))

  ##Hypothesis testing

  pos_int <- rep(NA, n_treat)
  names(pos_int) <- ad_int %>% filter(Treatment != "Control") %>% pull(Treatment)

  for (treat in 1:n_treat) {
    trteff_pos <- psample_int %>% pull(paste0("trteff[", treat, "]"))

    ## Extract the mean and SD of the treatment effect
    mean_trteff <- mean(trteff_pos)
    sd_trteff <- sd(trteff_pos)


    pos_int[treat] <- quantile(trteff_pos, (1 - gamma)) > 0
  }

  return(pos_int)

}

predict_success <- function(data,
                            N_planned,
                            gamma = 0.975,
                            num_chains, n.adapt, n.iter, perc_burnin) {

  n_arms <- length(unique(data$trt))
  n_treat <- n_arms - 1

  # Derive summary statistics at interim
  ad_int <- data  %>% group_by(Treatment) %>%
    summarize(n = n(),
              mean = mean(Y), var = var(Y),
              tau = n/var(Y))

  N <- rep(NA, n_arms)
  Treatments <- rep("", n_arms)
  Y <- data.frame(matrix(NA, ncol = n_arms, nrow = max(ad_int %>% pull(n))))
  colnames(Y) <- ad_int %>% pull(Treatment)

  for (k in seq(ncol(Y))) {
    N[k] <- (ad_int %>% pull(n))[k]
    Treatments[k] <- as.character((ad_int %>% pull(Treatment))[k])
    Y[1:N[k], k] <- data %>% filter(Treatment == Treatments[k]) %>% pull(Y)
  }


  jags_data <- list(n = N,
                    Y = Y,
                    mu = c(0, 0.4, 0.5),
                    tau = c(1/1, 1/1, 1/1),
                    n_pln = N_planned,
                    n_arms = n_arms)

  # Parameters to monitor
  params <- c("pred_trteff")

  # Analyse the interim data
  jags_model_int <- jags.model(file = file.path("./inst/jags/ppos_diff_means_multiarm.jags"),
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

  ##Hypothesis testing
  sign <- data.frame(matrix(NA, ncol = n_arms, nrow = max(ad_int %>% pull(n))))
  colnames(sign) <- Treatments
  colnames(sign)[1] <- "N_per_arm"

  for (treat in 2:n_arms) {
    for (ni in N_planned:max(ad_int %>% pull(n))) {
      trteff <- psample_int %>% pull(paste0("pred_trteff[", ni, ",", treat, "]"))
      sign[paste0("N_",ni), "N_per_arm"] <- ni
      sign[paste0("N_",ni), treat] <- ifelse(quantile(trteff, (1 - gamma)) > 0, 1, 0)
    }
  }

  # Delete rows with no data
  sign <- sign %>% filter(N_per_arm >= N_planned)

  return(sign)
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


calc_cp <- function(mu_c, mu_t, sd_c, sd_t, margin = 0,
                       n_int_c, n_pln_c, n_max_c,
                       kappa, alpha) {
  power_trial_start <- data.frame(matrix(NA, nrow = ntreat, ncol = 3))
  colnames(power_trial_start) <- c("N_int", "N_pln", "N_max")
  rownames(power_trial_start) <- names(mu_t)

  for (treat in seq(mu_t)) {
    power_trial_start[treat, "N_int"] <- cp_diff_means(n1 = n_int_c, # Sample size for the control arm
                                                        mu1 = mu_t[treat],
                                                        mu2 = mu_c,
                                                        margin = margin, #non-inferiority/superiority margin
                                                        sd1 = sd_t[treat],
                                                        sd2 = sd_c,
                                                        kappa = kappa, # allocation ratio for the active treatment
                                                        alpha = alpha,
                                                        alternative = "superiority")$power
    power_trial_start[treat, "N_pln"] <- cp_diff_means(n1 = n_pln_c, # Sample size for the control arm
                                                        mu1 = mu_t[treat],
                                                        mu2 = mu_c,
                                                        margin = margin, #non-inferiority/superiority margin
                                                        sd1 = sd_t[treat],
                                                        sd2 = sd_c,
                                                        kappa = kappa, # allocation ratio for the active treatment
                                                        alpha = alpha,
                                                        alternative = "superiority")$power
    power_trial_start[treat, "N_max"] <- cp_diff_means(n1 = n_max_c, # Sample size for the control arm
                                                        mu1 = mu_t[treat],
                                                        mu2 = mu_c,
                                                        margin = margin, #non-inferiority/superiority margin
                                                        sd1 = sd_t[treat],
                                                        sd2 = sd_c,
                                                        kappa = kappa, # allocation ratio for the active treatment
                                                        alpha = alpha,
                                                        alternative = "superiority")$power

  }
  power_trial_start
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
