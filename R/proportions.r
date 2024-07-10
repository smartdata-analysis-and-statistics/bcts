#' Sample size determination for Pearson Chi-square test
#'
#' @param p_treat Expected response probability in treated patients
#' @param p_control Expected response probability in control patients
#' @param delta Margin for deciding superiority (delta > 0) or non-inferiority (delta < 0)
#' @param kappa Treatment allocation ratio
#' @param alpha Desired level of significance
#' @param power Target power
#' @param test Type of test
#' @param design Study design
#'
#' @description
#' The problem of testing non-inferiority and superiority can be unified by the following hypotheses:
#' H0: p_treat - p_control <= delta
#' H1: p_treat - p_control > delta
#'
#' @source
#' Chow SC, Shao J Wang, Hansheng, Lokhnygina Y. Sample size calculations in clinical research. 2018.
#'
#' @examples
#' N_two_sample_proportions(p_treat = 0.8,
#'                          p_control = 0.8,
#'                          delta = -0.20,
#'                          kappa = 1,
#'                          alpha = 0.05,
#'                          power = 0.80)
#' @return A list
#' @export

N_two_sample_proportions <- function(p_treat = NA,
                                     p_control = NA,
                                     delta = 0,
                                     kappa = 1,
                                     alpha = 0.05,
                                     power = 0.8,
                                     test = "noninferiority/superiority",
                                     design = "parallel") {

  if (power < 0 | power > 1) {
    stop("Invalid value for 'power'")
  }
  if (alpha <= 0 | alpha >= 1) {
    stop("Invalid value for 'alpha'")
  }
  if (kappa <= 0) {
    stop("Invalid value for 'kappa'")
  }
  if (p_treat < 0 | p_treat > 1 | p_control < 0 | p_control > 1) {
    stop("The expected outcome proportions should be between 0 and 1")
  }

  if (test == "noninferiority/superiority" & design == "parallel") {
    z.alpha <-  stats::qnorm(alpha, mean = 0, sd = 1)
    z.beta <- stats::qnorm(1 - power, mean = 0, sd = 1)

    n_control <- (p_treat * (1 - p_treat)/kappa +
                            p_control * (1 - p_control)) * ((z.alpha + z.beta)/(p_treat - p_control - delta))^2
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
#' @param n_control Number of patients in the control group
#' @param p_treat Expected proportion in the treated group
#' @param p_control Expected proportion in the control group
#' @param delta Margin
#' @param tar Treatment allocation ratio
#' @param alpha Alpha
#' @param nsim Number of simulations to conduct
#'
#' @import stats
#'
#' @return List with result
#' @export
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
    pval <- pnorm(z_t, mean = 0, sd = 1, lower.tail = FALSE)

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

