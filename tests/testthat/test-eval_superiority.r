test_that("Test type-I error control Bayesian analysis", {
  nsim <- 10000
  alpha <- 0.05

  sig <- data.frame(treat1 = rep(NA, nsim), treat2 = NA)

  z_crit <- qnorm(1 - alpha)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 100, mean = c(0, 0, 0), sd = c(1, 1, 2),
                         trtnames = c("control", "treat1", "treat2"),
                         block.sizes = 1)

    result <- eval_superiority(data = ds, margin = 0, gamma = 1 - alpha,
                               method = "bayes", num_chains = 4,
                               n.adapt = 500, n.iter = 1000, perc_burnin = 0.2)

    # Ensure the Treatment column is treated as a character
    treatments <- as.character(result %>% pull(Treatment))

    sig[i, treatments] <- result$rejectH0
  }

  treat1.ci <- binom.confint(sum(sig$treat1), n = nsim, conf.level = 0.95, methods = "exact")
  treat2.ci <- binom.confint(sum(sig$treat2), n = nsim, conf.level = 0.95, methods = "exact")

  expect_true(alpha >= treat1.ci$lower & alpha <= treat1.ci$upper,
              info = paste0("Type-I error for treat1 is not maintained (", round(treat1.ci$lower,4), "; ", round(treat1.ci$upper,4)))
  expect_true(alpha >= treat2.ci$lower & alpha <= treat2.ci$upper,
              info = paste0("Type-I error for treat2 is not maintained (", round(treat2.ci$lower,4), "; ", round(treat2.ci$upper,4)))

})

test_that("Test type-I error control frequentist analysis", {
  nsim <- 10000
  alpha <- 0.05

  sig <- data.frame(treat1 = rep(NA, nsim), treat2 = NA)

  z_crit <- qnorm(1 - alpha)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 100, mean = c(0, 0, 0), sd = c(1, 1, 2),
                         trtnames = c("control", "treat1", "treat2"),
                         block.sizes = 1)

    result <- eval_superiority(data = ds, margin = 0, gamma = 1 - alpha, method = "mcmc")

    # Ensure the Treatment column is treated as a character
    treatments <- as.character(result %>% pull(Treatment))

    sig[i, treatments] <- result$rejectH0
  }

  treat1.ci <- binom.confint(sum(sig$treat1), n = nsim, conf.level = 0.95, methods = "exact")
  treat2.ci <- binom.confint(sum(sig$treat2), n = nsim, conf.level = 0.95, methods = "exact")

  expect_true(alpha >= treat1.ci$lower & alpha <= treat1.ci$upper,
              info = paste0("Type-I error for treat1 is not maintained (", round(treat1.ci$lower,4), "; ", round(treat1.ci$upper,4)))
  expect_true(alpha >= treat2.ci$lower & alpha <= treat2.ci$upper,
              info = paste0("Type-I error for treat2 is not maintained (", round(treat2.ci$lower,4), "; ", round(treat2.ci$upper,4)))

})

test_that("Test equivalence of Bayesian and frequentist approach", {

})
