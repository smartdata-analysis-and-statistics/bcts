test_that("Test type-I error control of test statistic", {
  require(binom)

  nsim <- 10000
  alpha <- 0.02

  sig = rep(NA, nsim)

  z_crit <- qnorm(1 - alpha)

  for (i in 1:nsim) {

    # Evaluate significance using frequentist method
    Yc <- rnorm(80, mean = 0, sd = 1)
    Yt <- rnorm(80, mean = 0, sd = 1)
    sepooled <- sqrt(var(Yc)/length(Yc) + var(Yt)/length(Yt))
    z_freq <- (mean(Yt) - mean(Yc))/sepooled
    sig[i] <- z_freq > z_crit
  }

  type1.ci <- binom.confint(sum(sig), n = nsim, conf.level = 0.95, methods = "exact")

  # Check if alpha is within the confidence interval
  expect_true(alpha >= type1.ci$lower & alpha <= type1.ci$upper,
              info = "The confidence interval does not contain alpha.")

})

