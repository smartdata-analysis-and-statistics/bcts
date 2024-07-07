library(dplyr)


source("R/simulation.r")

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
#' @param alternative
#'
#' @description
#' For superiority, H0 is defined as mu1 - mu2 <= margin.
#'
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
#' @param method Analysis method. Choose 'mcmc' or 'bayes'
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
#' @param quiet Display messaging output
#'
#' @description
#' The null hypothesis can be rejected when the posterior probability that $\\mu_t$ - $\\mu_c$ exceeds a high probability threshold,
#' i.e. $Pr(\\mu_t - \\mu_c > 0|D) > $\\gamma$)$.
#'
#'
#' @return An object of class "bcts"
#' @importFrom dplyr %>%
#' @importFrom dplyr slice_head mutate select filter add_row pull
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
bcts <- function(n_per_arm_int = 20,
                 n_per_arm_pln = 65,
                 n_per_arm_max = 80,
                 block.sizes = 1, # block size for randomization
                 mu_c = 0,
                 sd_c = 1,
                 mu_t = c("Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
                 sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1),
                 prioritize_low_trteff = TRUE,
                 gamma = 0.975, #0.98   %
                 th.fut = 0.2, ##
                 th.eff = 0.9, ##
                 th.prom = 0.5, ##
                 method = "mcmc",
                 nsim = 1000, # no. of simulatoins at trial start
                 num_chains = 4,
                 n.iter = 5000, #
                 n.adapt = 500,
                 perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                 progress.bar = "text") {


  # Sample the random seeds
  seeds <- sample(seq(nsim*10,nsim*20), size = nsim)

  # Select treatment with highest dose
  trt_high <- order(mu_t, decreasing = TRUE)[1] + 1
  trt_low <- order(mu_t, decreasing = FALSE)[1] + 1


  n_treat <- length(mu_t)
  n_arms <- n_treat + 1

  # Calculate the sample sizes
  n_pln <- 2*n_per_arm_pln + (length(mu_t) - 1)*n_per_arm_int
  n_max <- 2*n_per_arm_max + (length(mu_t) - 1)*n_per_arm_int

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

  results_interim <- results_final <- data.frame(sim = numeric(),
                           Treatment = character(),
                           trt = numeric(),
                           n = numeric(),
                           N = numeric(),
                           est_lower = numeric(),
                           est = numeric(),
                           est_upper = numeric(),
                           est_se = numeric(),
                           Z_test = numeric(),
                           Z_crit = numeric(),
                           I = numeric(),
                           rejectH0 = logical(),
                           ppos = numeric())
  names(results_interim)[names(results_interim) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(results_interim)[names(results_interim) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))
  names(results_final)[names(results_final) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(results_final)[names(results_final) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))

  if (progress.bar == "text") {
    pb <- utils::txtProgressBar(min = 0, max = nsim, initial = 0)
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


    # Evaluate PPOS at the planned sample size
    dat_pln <- dat_max %>%
      dplyr::slice_head(n = n_per_arm_pln*n_arms) %>%
      dplyr::mutate(Y = ifelse(id > n_per_arm_int*n_arms, NA, Y))

    result_pln <- eval_superiority(dat_pln, margin = 0, gamma = gamma,
                                   method = method, num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)

    results_interim <- results_interim %>% add_row(cbind(sim = i, result_pln))

    if (max(result_pln$ppos) < th.fut) {
      simresults$fut.trig[i] <- TRUE
      simresults$inc.ss[i] <- FALSE
      simresults$n_per_arm_fin[i] <- n_per_arm_pln
      simresults$sel.dose[i] <- result_pln$trt[order(result_pln$ppos, decreasing = TRUE)[1]]
      simresults$sel.dose.ppos[i] <- result_pln$ppos[order(result_pln$ppos, decreasing = TRUE)[1]]
    } else if (max(result_pln$ppos) >= th.prom & max(result_pln$ppos) < th.eff) {
      simresults$fut.trig[i] <- FALSE
      simresults$inc.ss[i] <- TRUE

      # Evaluate predictive power for each sample size
      ppos_by_n <- data.frame(n_per_arm_new = numeric(),
                              Treatment = character(),
                              trt = numeric(),
                              n = numeric(),
                              N = numeric(),
                              ppos = numeric())

      for (n_per_arm_new in (n_per_arm_pln + 1):n_per_arm_max) {

        dat_pln_new <- dat_max %>% slice_head(n = n_per_arm_new*n_arms) %>%
          mutate(Y = ifelse(id > n_per_arm_int*n_arms, NA, Y))

        result_pln_new <- eval_superiority(data = dat_pln_new, margin = 0,
                                           gamma = gamma, method = method,
                                           num_chains = num_chains,
                                           n.adapt = n.adapt, n.iter = n.iter,
                                           perc_burnin = perc_burnin)

        ppos_by_n <- ppos_by_n %>% add_row(cbind(n_per_arm_new, result_pln_new %>% select(Treatment, trt, n, N, ppos)))
      }

      # Select the lowest sample size resulting in a PPOS >= th.eff
      filtered_ppos <- subset(ppos_by_n, ppos >= th.eff)

      if (nrow(filtered_ppos) > 0) {
        sel_row <- filtered_ppos[which.min(filtered_ppos$n_per_arm_new),]
        simresults$sel.dose[i] <- sel_row$trt
        simresults$sel.dose.ppos[i] <- sel_row$ppos
        simresults$n_per_arm_fin[i] <- sel_row$n_per_arm_new
      } else {
        # Increasing the sample size did not help to attain target power
        # Keep the maximum sample size
        simresults$sel.dose[i] <- result_pln$trt[order(result_pln$ppos, decreasing = TRUE)[1]]
        simresults$sel.dose.ppos[i] <- result_pln$ppos[order(result_pln$ppos, decreasing = TRUE)[1]]
        simresults$n_per_arm_fin[i] <- n_per_arm_max
      }
    } else {
      simresults$fut.trig[i] <- FALSE
      simresults$inc.ss[i] <- FALSE
      simresults$n_per_arm_fin[i] <- n_per_arm_pln

      # If both predictive powers exceed 0.9 then
      # choose Low dose otherwise choose the dose with higher
      # predictive power
      if (sum((result_pln$ppos >= th.eff)) >= 2 & prioritize_low_trteff) {
        simresults$sel.dose[i] <- result_pln$trt[order(mu_t[which(result_pln$ppos >= th.eff)], decreasing = FALSE)[1]]
        simresults$sel.dose.ppos[i] <- result_pln$ppos[order(mu_t[which(result_pln$ppos >= th.eff)], decreasing = FALSE)[1]]
      } else {
        simresults$sel.dose[i] <- result_pln$trt[order(result_pln$ppos, decreasing = TRUE)[1]]
        simresults$sel.dose.ppos[i] <- result_pln$ppos[order(result_pln$ppos, decreasing = TRUE)[1]]
      }
    }

    # Keep the reference treatment and the selected treatment
    dat_fin <- dat_max %>%
      filter(trt %in% c(1,simresults$sel.dose[i])) %>%
      dplyr::slice_head(n = simresults$n_per_arm_fin[i]*2)

    result_fin <- eval_superiority(data = dat_fin, margin = 0, gamma = gamma,
                                   method = method, num_chains = num_chains,
                                   n.adapt = n.adapt, n.iter = n.iter,
                                   perc_burnin = perc_burnin)

    results_final <- results_final %>% add_row(cbind(sim = i, result_fin))

    simresults$sig_fin[i] <- result_fin %>% pull(rejectH0)

    if (progress.bar == "text") {
      utils::setTxtProgressBar(pb, i)
    }
  }

  if (progress.bar == "text") {
    close(pb)
  }

  out.n <- list(interim = n_per_arm_int*n_arms,
                final = mean(simresults$n_per_arm_fin*2) + n_per_arm_int*(n_arms - 2))

  out <- list(mu = c(mu_c, mu_t),
              sigma = c(sd_c, sd_t),
              n = out.n,
              gamma = gamma,
              power =  binom.confint(x = sum(simresults$sig_fin), n = nsim, methods = "exact"),
              fut.trig = binom.confint(x = sum(simresults$fut.trig), n = nsim, methods = "exact"),
              inc.ss = binom.confint(x = sum(simresults$inc.ss), n = nsim, methods = "exact"),
              sel.low = binom.confint(x = sum(simresults$sel.dose == trt_low), n = nsim, methods = "exact"),
              sel.high = binom.confint(x = sum(simresults$sel.dose == trt_high), n = nsim, methods = "exact"),
              simresults = simresults,
              results_interim = results_interim,
              results_final = results_final
              )

  class(out) <- "bcts"
  return(out)
}



#' Evaluate the posterior probability of treatment superiority
#'
#' @param data Trial data
#' @param margin Superiority margin (default: 0)
#' @param gamma Level to declare success (default: 0.975)
#' @param num_chains The number of parallel chains for the model
#' @param n.adapt The number of iterations for adaptation.
#' @param n.iter 	The number of iterations to monitor
#' @param perc_burnin The percentage of iterations to keep for burn-in
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
                             n.adapt = 500, n.iter = 1000, perc_burnin = 0.2) {

  if (method == "bayes") {
    result <- eval_superiority_bayes(data = data, margin = margin, gamma = gamma,
                                     num_chains = num_chains, n.adapt = n.adapt,
                                     n.iter = n.iter, perc_burnin = perc_burnin)
  } else if (method == "mcmc") {
    result <- eval_superiority_mcmc(data = data, margin = margin, gamma = gamma)
  } else {
    stop("Unknown method. Please use 'bayes' or 'mcmc'.")
  }

  return(result)
}


#' Bayesian method to evaluate superiority
#'
#' @param data
#' @param margin
#' @param gamma
#' @param num_chains
#' @param n.adapt
#' @param n.iter
#' @param perc_burnin
#'
#' @importFrom stats qnorm
#'
#' @return
#'
eval_superiority_bayes <- function(data, margin, gamma, num_chains, n.adapt, n.iter, perc_burnin) {

  treatments <- unique(data$trt)

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

  trt_est <- data.frame("Treatment" = character(),
                        "trt" = numeric(),
                        "n" = numeric(),
                        "N" = numeric(),
                        "est_lower" = numeric(),
                        "est" = numeric(),
                        "est_upper"  = numeric(),
                        "est_se" = numeric(), # Standard error of the mean
                        "Z_test" = numeric(),
                        "Z_crit" = numeric(),
                        "I" = numeric(), # Information level
                        "rejectH0" = logical(), # Reject H0?
                        "ppos" = numeric())  # Trial success based on observed data (i.e.,)

  for (treat in treatments[treatments != 1]) {
    mudiff <- psample_int %>% pull(paste0("trteff[", treat, "]"))
    PPOS <- psample_int %>% pull(paste0("ppos[", treat, "]"))
    trtc <-  (data %>% filter(trt == treat) %>% pull(Treatment))[1]
    n <-  nrow(data %>% filter(trt %in% c(1, treat) & !is.na(Y)))
    N <-  nrow(data %>% filter(trt %in% c(1, treat)))

    trt_est <- trt_est %>% add_row(data.frame("Treatment" = trtc,
                                              "trt" = as.numeric(treat),
                                              "n" = n,
                                              "N" = N,
                                              "est_lower" = stats::quantile(mudiff, 1 - gamma),
                                              "est" = mean(mudiff),
                                              "est_upper" = stats::quantile(mudiff, gamma),
                                              "est_se" = sd(mudiff),
                                              "Z_test" = mean(mudiff)/sd(mudiff),
                                              "Z_crit" = qnorm(gamma),
                                              "I" = 1/var(mudiff),
                                              "rejectH0" = stats::quantile(mudiff, (1 - gamma)) > 0,
                                              "ppos" = mean(PPOS)))
  }

  # Assign the dynamic column names
  names(trt_est)[names(trt_est) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(trt_est)[names(trt_est) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))

  return(trt_est)
}

eval_superiority_mcmc <- function(data, margin, gamma) {

  treatments <- unique(data$trt)

  trt_est <- data.frame("Treatment" = character(),
                        "trt" = numeric(),
                        "n" = numeric(),
                        "N" = numeric(),
                        "est_lower" = numeric(),
                        "est" = numeric(),
                        "est_upper"  = numeric(),
                        "est_se" = numeric(), # Standard error of the mean
                        "Z_test" = numeric(),
                        "Z_crit" = numeric(),
                        "I" = numeric(), # Information level
                        "rejectH0" = logical(), # Reject H0?
                        "ppos" = numeric())  # Trial success based on observed data (i.e.,)

  for (treat in treatments[treatments != 1]) {
    trtc <-  (data %>% filter(trt == treat) %>% pull(Treatment))[1]
    n <-  nrow(data %>% filter(trt %in% c(1, treat) & !is.na(Y)))
    N <-  nrow(data %>% filter(trt %in% c(1, treat)))
    Yc <- data %>% filter(trt == 1 & !is.na(Y)) %>% pull(Y)
    Yt <- data %>% filter(trt == treat & !is.na(Y)) %>% pull(Y)

    test_freq <- pow_diff_means(n1 = length(Yt), n2 = length(Yc), N = N,
                                mu1 = mean(Yt), mu2 = mean(Yc), margin = 0,
                                sd1 = sd(Yt), sd2  = sd(Yc),
                                alpha = 1 - gamma,
                                alternative = "superiority")

    trt_est <- trt_est %>% add_row(data.frame("Treatment" = trtc,
                                              "trt" = as.numeric(treat),
                                              "n" = n,
                                              "N" = N,
                                              "est_lower" = test_freq$diff_lower,
                                              "est" = test_freq$diff,
                                              "est_upper" = test_freq$diff_upper,
                                              "est_se" = test_freq$diff_se,
                                              "Z_test" = test_freq$Z_test,
                                              "Z_crit" = test_freq$Z_crit,
                                              "I" = test_freq$I,
                                              "rejectH0" = test_freq$rejectH0,
                                              "ppos" = test_freq$ppos))
  }

  # Assign the dynamic column names
  names(trt_est)[names(trt_est) == "est_lower"] <- paste0("est_", sub("0\\.", "", 1 - gamma))
  names(trt_est)[names(trt_est) == "est_upper"] <- paste0("est_", sub("0\\.", "", gamma))

  return(trt_est)
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



