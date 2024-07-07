test_that("Test type-I error control frequentist analysis", {
  nsim <- 10000
  alpha <- 0.05
  nC <- nT <- 100
  N <- 500
  sig <- rep(FALSE, nsim)

  for (i in 1:nsim) {

    # Evaluate significance using frequentist method
    Yc <- rnorm(nC, mean = 0, sd = 1)
    Yt <- rnorm(nT, mean = 0, sd = 1)

    test <- pow_diff_means(n1 = nT, n2 = nC, N = N, mu1 = mean(Yt),
                           mu2 = mean(Yc), margin = 0,
                           sd1 = sd(Yt), sd2 = sd(Yc),
                           alpha = alpha, alternative = "superiority")

    sig[i] <- test$rejectH0
  }

  type1.ci <- binom.confint(sum(sig), n = nsim, conf.level = 0.95, methods = "exact")

  # Check if alpha is within the confidence interval
  expect_true(alpha >= type1.ci$lower & alpha <= type1.ci$upper,
              info = paste0("The confidence interval does not contain alpha.", round(type1.ci$lower,4), "; ", round(type1.ci$upper,4)))

})
