library(dplyr)


#' Title
#'
#' @param n_0
#' @param mu
#' @param sd
#' @param trt
#' @param delta
#' @param tar
#' @param alpha
#' @param power  Choose between "marginal", "disjunctive" or "conjunctive"
#' @param correction
#' @param nsim
#'
#' @return
#' @export
#'
#' @examples
power_diff_means <- function(n_0, # Sample size for the control arm
                             mu = c(0, 0, 0), # expected effect sizes
                             sd = c(1, 1, 1), # expected SD
                             trt = NULL,
                             delta = 0, #non-inferiority/superiority margin
                             tar = NULL, # Vector with allocation ratios
                             alpha = 0.05,
                             power = "marginal",
                             correction = "bonferroni",
                             nsim = 10000) {

  n_0 <- ifelse(n_0 < 1, 1, ceiling(n_0))
  K <- length(mu) - 1 # Number of active arms

  if (is.null(correction)) {
    correction <- "none"
  }
  if (is.null(trt)) {
    # Specify treatment names
    trt <- c("control", paste0("treatment_", seq(K)))
  }
  if (is.null(tar)) {
    tar <- rep(1, (K + 1))
  }
  n <- n_0 * tar
  names(n) <- trt

  # Determine the significance level
  if (K == 1) {
    if (correction != "none") message("No correction for multiplicity is needed and thus will not be applied!")
    gamma <- alpha
  } else if (correction == "bonferroni") {
    gamma <- alpha/K
  } else if (correction == "Sidak") {
    gamma <- 1 - (1 - alpha)^(1/K)
  } else {
    if (correction == "none") warning("No correction for multiplicity was specified!")
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
  names(n_sign_per_comparison) <- trt[-1]

  # Estimate each power
  powers <- data.frame(matrix(NA, nrow = 1, ncol = K + 2))
  colnames(powers) <- c("disjunctive", "conjunctive", paste0("marginal_", trt[-1]))

  powers$disjunctive <- sum(n_sign > 0)/nsim
  powers$conjunctive = sum(n_sign == K)/nsim

  act_trt <- trt[-1]

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
  names(dat) <- paste0("z_",trt[-1])

  out <- list(power = est.power,
              powers = powers,
              n_0 = n_0,
              n = n,
              N = sum(n),
              K = K, # Number of active treatments
              alpha = alpha,
              power_type = power,
              p_crit = gamma,
              correction = correction,
              z_test_lower = qnorm(1 - gamma))
  return(out)
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
