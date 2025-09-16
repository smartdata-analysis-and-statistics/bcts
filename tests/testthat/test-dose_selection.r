test_that("Dose selection at interim (frequentist)", {
  nsim   <- 1000
  alpha  <- 0.05
  trt_ref <- "Placebo"
  trt_act <- c("Drug A", "Drug B")
  n_pln  <- 20 + 65*2
  th.fut <- 0.2
  th.eff <- 0.9

  rejectH0_freq <- matrix(FALSE, nrow = nsim, ncol = length(trt_act),
                          dimnames = list(NULL, trt_act))

  for (i in seq_len(nsim)) {
    dat_int <- sim_rct_normal(
      n    = 20*3,
      mean = c("Placebo" = 0, "Drug A" = 0, "Drug B" = 0),
      sd   = c("Placebo" = 1, "Drug A" = 1, "Drug B" = 1)
    )

    res <- dose_selection(
      dat_int    = dat_int,
      n_pln      = n_pln,
      trt_ref    = trt_ref,
      trt_active = trt_act,
      gamma      = 1 - alpha/2,
      th.fut     = th.fut,
      th.eff     = th.eff,
      method     = "mcmc"
    )

    rej <- res$rejectH0
    # normalize to named logical vector over trt_act
    if (is.logical(rej) && length(rej) == 1L) rej <- setNames(rep(rej, length(trt_act)), trt_act)
    if (is.null(names(rej))) names(rej) <- trt_act
    rejectH0_freq[i, ] <- as.logical(rej[trt_act])
  }

  for (arm in trt_act) {
    ci <- mc_error_proportion(x = sum(rejectH0_freq[, arm]),
                              n = nsim, level = 1 - alpha)
    expect_true(alpha >= ci$lower & alpha <= ci$upper,
                info = paste0("Type-I error not maintained (freq, ", arm, "): ",
                              round(ci$lower, 4), "–", round(ci$upper, 4)))
  }
})

test_that("Dose selection at interim (Bayesian)", {
  skip_if_not_installed("rjags")

  nsim   <- 1000
  alpha  <- 0.05
  trt_ref <- "Placebo"
  trt_act <- c("Drug A", "Drug B")
  n_pln  <- 20 + 65*2
  th.fut <- 0.2
  th.eff <- 0.9

  rejectH0_bayes <- matrix(FALSE, nrow = nsim, ncol = length(trt_act),
                           dimnames = list(NULL, trt_act))

  for (i in seq_len(nsim)) {
    dat_int <- sim_rct_normal(
      n    = 20*3,
      mean = c("Placebo" = 0, "Drug A" = 0, "Drug B" = 0),
      sd   = c("Placebo" = 1, "Drug A" = 1, "Drug B" = 1)
    )

    res <- dose_selection(
      dat_int    = dat_int,
      n_pln      = n_pln,
      trt_ref    = trt_ref,
      trt_active = trt_act,
      gamma      = 1 - alpha/2,
      th.fut     = th.fut,
      th.eff     = th.eff,
      method     = "bayes"
    )

    rej <- res$rejectH0
    if (is.logical(rej) && length(rej) == 1L) rej <- setNames(rep(rej, length(trt_act)), trt_act)
    if (is.null(names(rej))) names(rej) <- trt_act
    rejectH0_bayes[i, ] <- as.logical(rej[trt_act])
  }

  for (arm in trt_act) {
    ci <- mc_error_proportion(x = sum(rejectH0_bayes[, arm]),
                              n = nsim, level = 1 - alpha)
    expect_true(alpha >= ci$lower & alpha <= ci$upper,
                info = paste0("Type-I error not maintained (Bayes, ", arm, "): ",
                              round(ci$lower, 4), "–", round(ci$upper, 4)))
  }
})

