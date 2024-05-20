#' Sample size determination for Pearson Chi-square test
#'
#' @param alpha Desired level of significance
#' @param power Target power
#' @param p_treat Expected response probability in treated patients
#' @param p_control Expected response probability in control patients
#' @param delta Margin for deciding superiority (delta > 0) or non-inferiority (delta < 0)
#' @param kappa Treatment allocation ratio
#'
#' @description
#' The problem of testing non-inferiority and superiority can be unified by the following hypotheses:
#' H0: p_treat - p_control <= delta
#' H1: p_treat - p_control > delta
#'
#' @source
#' Chow SC, Shao J Wang, Hansheng, Lokhnygina Y. Sample size calculations in clinical research. 2018.
#'
#' @example
#' N_two_sample_proportions(p_treat = 0.8,
#'                          p_control = 0.8,
#'                          delta = -0.20,
#'                          kappa = 1,
#'                          alpha = 0.05,
#'                          power = 0.80)
#' @return
#' @export
#'
#' @examples
N_two_sample_proportions <- function(p_treat = NA,
                                     p_control = NA,
                                     delta = 0,
                                     kappa = 1,
                                     alpha = 0.05,
                                     power = 0.8,
                                     test = "noninferiority/superiority",
                                     design = "parallel") {

  if (test == "noninferiority/superiority" & design == "parallel") {
    z.alpha <-  qnorm(alpha, mean = 0, sd = 1)
    z.beta <- qnorm(1 - power, mean = 0, sd = 1)

    n_control <- ceiling((p_treat * (1 - p_treat)/kappa +
                            p_control * (1 - p_control)) * ((z.alpha + z.beta)/(p_treat - p_control - delta))^2)
    n_treat <- kappa * n_control

    return(list(
      H0 = paste0("p_treat - p_control <= ", delta),
      H1 =  paste0("p_treat - p_control > ", delta),
      power = power,
      alpha = alpha,
      N = data.frame(n.control = n_control, n.treat = n_treat, n.total = n_control + n_treat)))
  }

}

#' Calculate the power for Fisher Test
#'
#' @param n_control
#' @param p_treat
#' @param p_control
#' @param delta
#' @param tar
#' @param alpha
#' @param nsim Number of simulations to conduct
#'
#' @return
#' @export
#'
#' @examples
power_fisher <- function(n_control,
                         p_treat = NA,
                         p_control = NA,
                         delta,
                         tar = 1,
                         alpha = 0.05,
                         nsim = 10000) {

  n_control <- ifelse(n_control < 1, 1, ceiling(n_control))
  n_treat <- n_control * tar

  sign <- rep(NA, nsim)

  for (i in seq(nsim)) {
    y_c = rbinom(n = n_control, prob = p_control, size = 1)
    y_t = rbinom(n = n_treat, prob = p_treat, size = 1)

    obs_p_t <- mean(y_t)
    obs_p_c <- mean(y_c)

    var_p_t <- obs_p_t*(1 - obs_p_t)/length(y_t)
    var_p_c <- obs_p_c*(1 - obs_p_c)/length(y_c)

    # Derive the z-score
    z_t <- (obs_p_t - obs_p_c - delta)/sqrt(var_p_t + var_p_c)

    # The p-value is determined as the probability of observing data as
    # extreme as or more extreme than the data observed, assuming that the
    # null hypothesis is true.
    # Pr(obs_p_t - obs_p_c - delta > 0 | p_treat - p_control <= delta)
    # Pr(z > 0 | p_treat - p_control <= delta)
    pval <- pnorm(z, mean = 0, sd = 1, lower.tail = FALSE)

    # Determine if the p-value is statistically significant
    sign[i] <- z_t > qnorm(1 - alpha) # pval < alpha
  }

  n <- data.frame(n_control = n_control,
                  n_treat = n_treat,
                  n_total = n_control + n_treat,
                  power = mean(sign))

  return(list(power = n$power,
              n = n))
}


#' Calculate the power for detecting
#'
#' @param n_control
#' @param mu_treat
#' @param mu_control
#' @param sd_treat
#' @param sd_control
#' @param delta
#' @param tar
#' @param alpha
#' @param nsim
#'
#' @return
#' @export
#'
#' @examples
power_diff_means <- function(n_control,
                             mu_treat = NA,
                             mu_control = NA,
                             sd_treat = NA,
                             sd_control = NA,
                             delta = 0,
                             tar = 1,
                             alpha = 0.05,
                             nsim = 10000) {

  n_control <- ifelse(n_control < 1, 1, ceiling(n_control))
  n_treat <- n_control * tar

  sign <- rep(NA, nsim)

  for (i in seq(nsim)) {
    y_c <- rnorm(n = n_control, mean = mu_control, sd = sd_control)
    y_t <- rnorm(n = n_treat, mean = mu_treat, sd = sd_treat)

    # Observed response probabilities
    m_t <- mean(y_t)
    m_c <- mean(y_c)

    # Error variance of the observed response probabilities
    var_pooled <- (1/(n_control + n_treat - 2)) * (sum((y_c - m_c)**2) + sum((y_t - m_t)**2))

    # Derive the z-score
    z <- (m_t - m_c - delta)/(sqrt(var_pooled)*sqrt(1/n_control + 1/n_treat))

    # Determine whether to reject H0
    sign[i] <- z > qnorm(1 - alpha)
  }

  n <- data.frame(n_control = n_control,
                  n_treat = n_treat,
                  n_total = n_control + n_treat,
                  power = mean(sign))

  return(list(power = n$power,
              n = n))
}

#' Calculate the power for noninferiority/superiority
#'
#' @param n_control
#' @param mu_treat
#' @param mu_control
#' @param sd_treat
#' @param sd_control
#' @param delta
#' @param tar
#' @param alpha
#' @param nsim
#'
#' @return
#' @export
#'
#' @examples
power_adapt_diff_means <- function(n_control_stage1,
                                   mu_treat = 0.5,
                                   mu_control = 0,
                                   sd_treat = 1,
                                   sd_control = 1,
                                   delta = 0,
                                   tar = 1,
                                   K = 2, # Number of stages
                                   r = c(1, 1.5),   # Allocation ratio stage 1, stage 2, ... stage K
                                   alpha = 0.05,
                                   correction = c("BF"), # BF = bonferroni
                                   nsim = 10000) {

  if (tar != 1) {
    stop("Alpha spending for unequal treatment groups is not implemented!")
  }
  if (correction == "BF") {
    alpha_c <- alpha/K
    C_p <-  qnorm(1 - alpha_c)
  } else if (correction == "PK") {
    require(ldbounds)
    C_p <- commonbounds(looks = K, iuse = "PK", alpha = alpha, sides = 1)
  }

  n_control_stage1 <- ifelse(n_control_stage1 < 1, 1, ceiling(n_control_stage1))
  n_treat_stage1 <- n_control_stage1 * tar

  n_control <- ceiling(n_control_stage1 * sum(r))
  n_treat <- ceiling(n_treat_stage1 * sum(r))

  results <- data.frame(n_control = rep(NA, nsim),
                        n_treat = rep(NA, nsim),
                        rejectH0 = rep(NA, nsim))

  for (i in seq(nsim)) {

    # Sample all data
    y_c <- rnorm(n = n_control, mean = mu_control, sd = sd_control)
    y_t <- rnorm(n = n_treat, mean = mu_treat, sd = sd_treat)

    for (k in seq(K)) {
      n_c_k <- ceiling(n_control_stage1 * sum(r[1:k]))
      n_t_k <- ceiling(n_treat_stage1 * sum(r[1:k]))
      n_k <- n_c_k + n_t_k

      y_c_k <- y_c[1:n_c_k]
      y_t_k <- y_t[1:n_t_k]

      # Observed response probabilities
      m_c_k <- mean(y_c_k)
      m_t_k <- mean(y_t_k)

      # Derive the z-score
      var_pooled <- (1/(n_c_k + n_t_k - 2)) * (sum((y_c_k - m_c_k)**2) + sum((y_t_k - m_t_k)**2))
      z_k <- (m_t_k - m_c_k - delta)/(sqrt(var_pooled)*sqrt(1/n_c_k + 1/n_t_k))

      #z_k <- (m_t_k - m_c_k - delta)/(sqrt(n_k*(var(y_c_k) + var(y_t_k))))

      # Determine whether to reject H0
      if (z_k > C_p & k < K) {
        # Stop, reject H0
        results$n_control[i] <- n_c_k
        results$n_treat[i] <- n_t_k
        results$rejectH0[i] <- 1
        break
      } else if (z_k > C_p & k == K) {
        # Stop, reject H0
        results$n_control[i] <- n_c_k
        results$n_treat[i] <- n_t_k
        results$rejectH0[i] <- 1
      } else if (z_k <= C_p & k == K) {
        # Accept H0
        results$n_control[i] <- n_c_k
        results$n_treat[i] <- n_t_k
        results$rejectH0[i] <- 0
      }
    }
  }


  n <- data.frame(n_control = n_control,
                  n_control_mean = mean(results$n_control),
                  n_treat = n_treat,
                  n_treat_mean = mean(results$n_treat),
                  n_total = n_control + n_treat,
                  n_total_mean = mean(results$n_control+results$n_treat),
                  power = mean(results$rejectH0))

  return(list(power = n$power,
              n = n))
}


power_fisher_pocock <- function(n_control,
                                p_treat = NA,
                                p_control = NA,
                                delta,
                                tar = 1,
                                alpha = 0.05,
                                r_hist,
                                n_hist,
                                var_err_hist = 0,
                                nsim = 10000) {

  n_control <- ifelse(n_control < 1, 1, ceiling(n_control))
  n_treat <- n_control * tar

  sign <- rep(NA, nsim)

  for (i in seq(nsim)) {
    y_r = rbinom(n = n_control, prob = p_control, size = 1)
    y_t = rbinom(n = n_treat, prob = p_treat, size = 1)

    obs_p_t <- mean(y_t)
    obs_p_r <- mean(y_r)
    obs_p_h <- r_hist/n_hist

    var_p_t <- obs_p_t*(1 - obs_p_t)/length(y_t)
    var_p_r <- obs_p_r*(1 - obs_p_r)/length(y_r)
    var_p_h <- obs_p_h*(1 - obs_p_h)/n_hist

    # Derive a pooled estimate for the control response
    mu_c <- (((var_p_h + var_err_hist) * obs_p_r) + (var_p_r * obs_p_h)) /
      (var_p_r + var_p_h + var_err_hist)

    var_c <- 1/((1/var_p_r) + (1/(var_err_hist + var_p_h)))

    # Derive the z-score
    z <- (obs_p_t - mu_c - delta)/sqrt(var_p_t + var_c)

    # The p-value is determined as the probability of observing data as
    # extreme as or more extreme than the data observed, assuming that the
    # null hypothesis is true.
    # Pr(obs_p_t - obs_p_c - delta > 0 | p_treat - p_control <= delta)
    # Pr(z > 0 | p_treat - p_control <= delta)
    pval <- pnorm(z, mean = 0, sd = 1, lower.tail = FALSE)

    # Determine if the p-value is statistically significant
    sign[i] <- pval < alpha
  }

  n <- data.frame(n_control = n_control,
                  n_treat = n_treat,
                  n_total = n_control + n_treat,
                  power = mean(sign))

  return(list(power = n$power,
              n = n))
}
