test_that("Test type-I error control", {
  require(binom)

  nsim <- 10000
  alpha <- 0.02

  sig = rep(NA, nsim)

  z_crit <- qnorm(1 - alpha)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 100, mean = c(0, 0), sd = c(1, 1),
                         trtnames = c("control", "treat"),
                         block.sizes = 1)

    # Evaluate significance using frequentist method
    Yc <- ds %>% filter(Treatment == "control") %>% pull(Y)
    Yt <- ds %>% filter(Treatment == "treat") %>% pull(Y)
    sepooled <- sqrt(var(Yc)/length(Yc) + var(Yt)/length(Yt))
    z_freq <- (mean(Yt) - mean(Yc))/sepooled
    sig[i] <- z_freq > z_crit
  }

  type1.ci <- binom.confint(sum(sig), n = nsim, conf.level = 0.95, methods = "exact")

  # Check if alpha is within the confidence interval
  expect_true(alpha >= type1.ci$lower & alpha <= type1.ci$upper,
              info = "The confidence interval does not contain alpha.")

})


test_that("Test group means", {

  nsim <- 1000

  mu_est <- data.frame(control = rep(NA, nsim), treat1 = NA, treat2 = NA)

  for (i in 1:nsim) {
    # generate trial data
    ds <- sim_rct_normal(n = 500, mean = c(0, 0.4, 0.5), sd = c(1, 1, 1),
                         trtnames = c("control", "treat1", "treat2"),
                         block.sizes = 1)

    means <- ds %>% group_by(Treatment) %>% summarize(mu = mean(Y))

    # Ensure the Treatment column is treated as a character
    treatments <- as.character(means %>% pull(Treatment))

    mu_est[i, treatments] <- means %>% pull(mu)
  }

  se <- apply(mu_est, 2, sd)/sqrt(nsim)


  # Check if alpha is within the confidence interval
  expect_true(alpha >= type1.ci$lower & alpha <= type1.ci$upper,
              info = "The confidence interval does not contain alpha.")

})
