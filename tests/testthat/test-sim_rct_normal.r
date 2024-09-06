require(binom)

test_that("Test type-I error control", {
  nsim <- 10000
  alpha <- 0.02

  sig = rep(NA, nsim)

  z_crit <- qnorm(1 - alpha)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 100, mean = c(0, 0), sd = c(1, 1),
                         trtnames = c("control", "treat"))

    # Evaluate significance using frequentist method
    Yc <- ds %>% dplyr::filter(Treatment == "control") %>% pull(Y)
    Yt <- ds %>% dplyr::filter(Treatment == "treat") %>% pull(Y)
    sepooled <- sqrt(var(Yc)/length(Yc) + var(Yt)/length(Yt))
    z_freq <- (mean(Yt) - mean(Yc))/sepooled
    sig[i] <- z_freq > z_crit
  }

  type1.ci <- binom.confint(sum(sig), n = nsim, conf.level = 0.95, methods = "exact")

  # Check if alpha is within the confidence interval
  expect_true(alpha >= type1.ci$lower & alpha <= type1.ci$upper,
              info = "The confidence interval does not contain alpha.")

})


test_that("Test distribution of simulated data", {

  nsim <- 1000

  mu_est <- data.frame(control = rep(NA, nsim), treat1 = NA, treat2 = NA)
  sigma_est <- data.frame(control = rep(NA, nsim), treat1 = NA, treat2 = NA)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 1000, mean = c(0, 0.4, 0.5), sd = c(1, 1, 2),
                         trtnames = c("control", "treat1", "treat2"))

    pop_param <- ds %>% group_by(Treatment) %>% summarize(mu = mean(Y), sigma = sd(Y))

    # Ensure the Treatment column is treated as a character
    treatments <- as.character(pop_param %>% pull(Treatment))

    mu_est[i, treatments] <- pop_param %>% pull(mu)
    sigma_est[i, treatments] <- pop_param %>% pull(sigma)
  }

  mu_qnt <- apply(mu_est, 2, quantile, c(0.05, 0.95))
  sigma_qnt <- apply(sigma_est, 2, quantile, c(0.05, 0.95))


  expect_true(0 >= mu_qnt["5%", "control"] & 0 <= mu_qnt["95%", "control"],
              info = "The mean outcome in treatment group 'control' is not controlled")
  expect_true(0.4 >= mu_qnt["5%", "treat1"] & 0.4 <= mu_qnt["95%", "treat1"],
              info = "The mean outcome in treatment group 1 is not controlled.")
  expect_true(0.5 >= mu_qnt["5%", "treat2"] & 0.5 <= mu_qnt["95%", "treat2"],
              info = "The mean outcome in treatment group 2 is not controlled")

  expect_true(1 >= sigma_qnt["5%", "control"] & 1 <= sigma_qnt["95%", "control"],
              info = "The mean outcome in treatment group 'control' is not controlled")
  expect_true(1 >= sigma_qnt["5%", "treat1"] & 1 <= sigma_qnt["95%", "treat1"],
              info = "The mean outcome in treatment group 1 is not controlled.")
  expect_true(2 >= sigma_qnt["5%", "treat2"] & 2 <= sigma_qnt["95%", "treat2"],
              info = "The mean outcome in treatment group 2 is not controlled")

})
