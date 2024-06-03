library(dplyr)


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
power_diff_means <- function(n_0, # Sample size for the control arm
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
  uniroot(function(x) (power_diff_means(n_0 = x,
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

bpower_diff_means <- function(n_0, #planned initial sample size for the control group
                              p_interim = 0.3, # The information fraction at which the interim analysis is carried out
                              mu_c = 0,
                              sd_c = 1,
                              mu_t = c(0.4, 0.5), # expected effect sizes
                              sd_t = c(1, 1), # expected SD
                              gamma = 0.975, # Success threshold
                              nsim = 1000,
                              num_chains = 4,
                              n.iter = 20000, # Number of MCMC iterations
                              n.adapt = 500,
                              perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                              quiet = TRUE) {
  require(rjags)

  ntreat <- length(mu_t)

  ## Calculate the interim sample size
  n_0_interim <- ceiling(n_0 * p_interim)


  # Parameters to monitor
  params <- c("mu_t", "mu_c", "trteff")

  ##############################################################################
  # Simulation results to monitor
  ##############################################################################

  # keep track of predictive power for each dose at interim
  pred_power_int <- matrix(nrow = nsim, ncol = ntreat)

  # Keep track of the best doses at interim
  best_dose <- rep(NA, nsim)

  # Keep track of hypothesis test at final stage
  hypothesis_testing <- rep(NA, nsim)

  ##############################################################################
  # Run simulation
  ##############################################################################
  for (i in 1:nsim) {

    # Simulate new data for the control group
    control_data <- rnorm(n_0_interim, mu_c, sd_c)

    # Simulate new data for the treatment groups
    treatment_data <- matrix(nrow = n_0_interim, ncol = ntreat)
    mean_t <- tau_t <- rep(NA, ntreat)

    for (treat in 1:ntreat) {
      treatment_data[,treat] <- rnorm(n_0_interim, mu_t[treat], sd_t[treat])
      mean_t[treat] <- mean(treatment_data[,treat])
      tau_t[treat] <- n_0_interim/var(treatment_data[,treat]) # 1/(SE(mean))**2
    }

    jags_data <- list(mean_c = mean(control_data),
                      tau_c = n_0_interim/var(control_data),
                      mean_t = mean_t,
                      tau_t = tau_t)

    # Analyse the interim data
    jags_model_int <- jags.model(file = file.path("./inst/jags/multiarm.jags"),
                                 data = jags_data,
                                 n.chains = num_chains,
                                 n.adapt = n.adapt,
                                 quiet = quiet)

    ##Burnin stage
    update(jags_model_int, n.iter = ceiling(n.iter*perc_burnin), progress.bar = "none")

    ##Sampling after burnin
    fit_jags_int <- coda.samples(jags_model_int, params,
                                 n.iter = floor(n.iter*(1 - perc_burnin)),
                                 progress.bar = "none")

    ## Derive posterior distributions at interim
    psample_interim <- as.data.frame(do.call(rbind, fit_jags_int))

    ## Calculate predictive powers each dose
    for (treat in 1:ntreat) {
      pred_power_int[i, treat] <- mean(psample_interim[[paste0("trteff[", treat, "]")]] > 0)
    }

    # Select dose with highest predictive power
    best_dose[i] <- order(pred_power_int[i,], decreasing = TRUE)[1]

    n_0_fin <- n_0

    ## Simulate the remaining control patients
    control_data_fin <- c(control_data, rnorm(n_0_fin - n_0_interim, mu_c, sd_c))

    # Simulate the remaining treated patients for the selected dose
    treatment_data_fin <- c(treatment_data[,best_dose[i]],  rnorm(n_0_fin - n_0_interim, mu_t[best_dose[i]], sd_t[best_dose[i]]))

    jags_data_fin <- list(mean_c = mean(control_data_fin),
                          tau_c = n_0_fin/var(control_data_fin),
                          mean_t = mean(treatment_data_fin),
                          tau_t = n_0_fin/var(treatment_data_fin))

    ## Specify final JAGS model
    jags_model_fin <- jags.model(file = file.path("./inst/jags/multiarm.jags"),
                                 data = jags_data_fin,
                                 n.chains = num_chains,
                                 n.adapt = n.adapt,
                                 quiet = quiet)

    ##Burnin stage
    update(jags_model_fin, n.iter = ceiling(n.iter*perc_burnin), progress.bar = "none")

    ##Sampling after burnin
    fit_jags_fin <- coda.samples(jags_model_fin, params,
                                 n.iter = floor(n.iter*(1 - perc_burnin)),
                                 progress.bar = "none")

    ## Derive posterior distributions at interim
    psample_fin <- as.data.frame(do.call(rbind, fit_jags_fin))

    ##Hypothesis testing
    hypothesis_testing[i] <- (quantile(psample_fin$trteff, 0.025) > 0)
  }





}

