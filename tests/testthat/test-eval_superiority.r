test_that("Test type-I error control for traditional RCT", {
  skip_if_not_installed("rjags")

  nsim <- 10000
  alpha <- 0.05

  sigf <- sigb <- rep(NA, nsim)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 100, mean = c(0, 0), sd = c(1, 1),
                         trtnames = c("control", "treat"))

    sigf[i] <- eval_superiority(data = ds, margin = 0, gamma = 1 - alpha,
                               method = "mcmc", trt_ref = "control")$rejectH0
    sigb[i] <- eval_superiority(data = ds, margin = 0, gamma = 1 - alpha,
                                method = "bayes", num_chains = 4,
                                n.adapt = 500, n.iter = 1000, perc_burnin = 0.2,
                                trt_ref = "control")$rejectH0
  }

  treat.ci <- mc_error_proportion(x = sum(sigf), n = nsim, level = 1 - alpha)
  expect_true(alpha >= treat.ci$lower & alpha <= treat.ci$upper,
              info = paste0("Type-I error for treat is not maintained for frequentist analysis (",
                            round(treat.ci$lower,4), "; ", round(treat.ci$upper,4), ")"))

  treat.ci <- mc_error_proportion(x = sum(sigb), n = nsim, level = 1 - alpha)
  expect_true(alpha >= treat.ci$lower & alpha <= treat.ci$upper,
              info = paste0("Type-I error for treat is not maintained for bayesian analysis (",
                            round(treat.ci$lower,4), "; ", round(treat.ci$upper,4), ")"))

  # Estimation method should not be related to statistical significance
  ds.test <- data.frame(method = 0, sign = as.numeric(sigf)) %>%
    add_row(method = 1, sign = as.numeric(sigb))

  fit <- glm(sign ~ method, data = ds.test, family = binomial)

  p_value <- summary(fit)$coefficients["method", "Pr(>|z|)"]

  expect_true(p_value > alpha,
              info = paste0("There is a relation between estimation method and statistical significance (slope = ",
                            round(summary(fit)$coefficients["method", "Estimate"],4), "; ",
                            round(summary(fit)$coefficients["method", "Std. Error"],4), ")"))
})

test_that("Test power for traditional RCT", {
  skip_if_not_installed("rjags")

  nsim <- 10000
  alpha <- 0.05

  sigf <- sigb <- rep(NA, nsim)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 100, mean = c(0, 0.5), sd = c(1, 1),
                         trtnames = c("control", "treat"))

    sigf[i] <- eval_superiority(data = ds, margin = 0, gamma = 1 - alpha,
                                method = "mcmc", trt_ref = "control")$rejectH0
    sigb[i] <- eval_superiority(data = ds, margin = 0, gamma = 1 - alpha,
                                method = "bayes", num_chains = 4,
                                n.adapt = 500, n.iter = 1000, perc_burnin = 0.2,
                                trt_ref = "control")$rejectH0
  }

  # Estimation method should not be related to statistical significance
  ds.test <- data.frame(method = 0, sign = as.numeric(sigf)) %>%
    add_row(method = 1, sign = as.numeric(sigb))

  fit <- glm(sign ~ method, data = ds.test, family = binomial)

  p_value <- summary(fit)$coefficients["method", "Pr(>|z|)"]

  expect_true(p_value > alpha,
              info = paste0("There is a relation between estimation method and statistical significance (slope = ",
                            round(summary(fit)$coefficients["method", "Estimate"],4), "; ",
                            round(summary(fit)$coefficients["method", "Std. Error"],4), ")"))
})


