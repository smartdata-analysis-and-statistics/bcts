require(rjags)

test_that("Test dose selection at interim", {
  nsim <- 1000
  alpha <- 0.05

  trt_ref <- "Placebo"
  n_pln <- 20 + 65*2
  th.fut <- 0.2
  th.eff <- 0.9
  th.prom <- 0.5


  rejectH0  <- data.frame(freq = logical(), bayes = logical())

  for (i in 1:nsim) {

    # Simulate trial data at interim stage
    dat_int <- sim_rct_normal(n = 20*3,
                              mean = c("Placebo" = 0, "Drug A" = 0, "Drug B" = 0),
                              sd = c("Placebo" = 1, "Drug A" = 1, "Drug B" = 1))

    resultf <- dose_selection(dat_int = dat_int,
                              n_pln = n_pln,
                              trt_ref = "Placebo",
                              trt_active = c("Drug A", "Drug B"),
                              gamma = 1 - alpha/2,
                              th.fut = th.fut, th.eff = th.eff,
                              method = "mcmc")

    resultb <- dose_selection(dat_int = dat_int,
                              n_pln = n_pln,
                              trt_ref = "Placebo",
                              trt_active = c("Drug A", "Drug B"),
                              gamma = 1 - alpha/2,
                              th.fut = th.fut, th.eff = th.eff,
                              method = "bayes")

    rejectH0 <- rejectH0 %>% add_row(data.frame(freq = resultf$rejectH0, bayes = resultb$rejectH0))
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
