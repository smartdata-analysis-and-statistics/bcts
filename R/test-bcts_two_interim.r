require(rjags)

test_that("Test RCT with two interim analyses", {
  nsim <- 10000
  alpha <- 0.05


  for (i in 1:nsim) {

    sim <- bcts_two_interim(n_int1 = 60, n_int2 = 150, n_max = 180,
                            mu =  c("Placebo" = 0, "Drug A" = 0, "Drug B" = 0),
                            sigma = c("Placebo" = 1, "Drug A" = 1, "Drug B" = 1),
                            trt_ref = "Placebo",
                            gamma = 0.975,
                            th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                            method = "mcmc", nsim = 1000)
  }

  # Test whether type-I error is controlled
  treat.ci <- mc_error_proportion(x = sum(rejectH0$freq), n = nsim, level = 1 - alpha)
  expect_true(alpha >= treat.ci$lower & alpha <= treat.ci$upper,
              info = paste0("Type-I error for treat is not maintained for frequentist analysis (",
                            round(treat.ci$lower,4), "; ", round(treat.ci$upper,4), ")"))

  treat.ci <- mc_error_proportion(x = sum(rejectH0$bayes), n = nsim, level = 1 - alpha)
  expect_true(alpha >= treat.ci$lower & alpha <= treat.ci$upper,
              info = paste0("Type-I error for treat is not maintained for Bayesian analysis (",
                            round(treat.ci$lower,4), "; ", round(treat.ci$upper,4), ")"))


})
