require(rjags)

test_that("Test type-I error control for traditional RCT", {
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



test_that("Test adaptive RCT", {
  nsim <- 10000
  alpha <- 0.05

  mu <- c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5)
  sigma <- c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1)
  trt_rank <- c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3)
  trt_ref <- "Placebo"

  # Extract the treatment names
  trt_names <- names(mu)
  names(mu) <- names(sigma) <- trt_names

  # Active treatment names
  trt_active <- setdiff(trt_names, trt_ref)

  sigf <- sigb <- rep(NA, nsim)

  inc.ss <- fut.trig  <- data.frame(freq = logical(), bayes = logical())
  n.final <- data.frame(freq = numeric(), bayes = numeric())

  for (i in 1:nsim) {

    # Simulate trial data at interim stage
    dat_int <- sim_rct_normal(n = 20*3,
                              mean = mu,
                              sd = sigma,
                              trtnames = trt_names)

    resultf <- eval_adaptive_trial(dat_int = dat_int,
                                  n_pln = 20 + 65*2,
                                  n_max = 20 + 80*2,
                                  mu = mu, sigma = sigma, trt_ref = trt_ref, trt_active = trt_active,
                                  trt_rank = trt_rank,
                                  prioritize_low_rank = TRUE, gamma = 0.98,
                                  th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                                  method = "mcmc")

    resultb <- eval_adaptive_trial(dat_int = dat_int,
                                   n_pln = 20 + 65*2,
                                   n_max = 20 + 80*2,
                                   mu = mu, sigma = sigma, trt_ref = trt_ref, trt_active = trt_active,
                                   trt_rank = trt_rank,
                                   prioritize_low_rank = TRUE, gamma = 0.98,
                                   th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                                   method = "bayes", num_chains = 4,
                                   n.iter = 5000, n.adapt = 500, perc_burnin = 0.2)

    inc.ss <- inc.ss %>% add_row(data.frame(freq = resultf$result$inc.ss, bayes = resultb$result$inc.ss))
    fut.trig <- fut.trig %>% add_row(data.frame(freq = resultf$result$fut.trig, bayes = resultb$result$fut.trig))
    n.final <- n.final %>% add_row(data.frame(freq = resultf$result$n.final, bayes = resultb$result$n.final))
  }

  test.fut.trig <- prop.test(x = colSums(fut.trig), n = rep(nsim, 2), alternative = "two.sided", correct = TRUE)
  expect_true(test.fut.trig$p.value > alpha,
              info = paste0("Futility triggering is not equally likely across the estimation methods (",
                            paste(colMeans(fut.trig), collapse = "; "), ")"))

  test.inc.ss <- prop.test(x = colSums(inc.ss), n = rep(nsim, 2), alternative = "two.sided", correct = TRUE)
  expect_true(test.inc.ss$p.value > alpha,
              info = paste0("Sample size increase is not equally likely across the estimation methods (",
                            paste(colMeans(inc.ss), collapse = "; "), ")"))

  test.n.final <- t.test(n.final$freq, n.final$bayes, var.equal = FALSE)
  expect_true(test.n.final$p.value > alpha,
              info = paste0("Final sample size is not equal across the estimation methods (",
                            paste(colMeans(n.final), collapse = "; "), ")"))

})
