test_that("C++ and R implementations yield similar power estimates", {
  skip_on_cran()

  set.seed(42)
  B <- 1000
  n_draws <- 5000

  # Common settings
  settings <- list(
    B = B,
    p_c = 0.85,
    p_t = 0.85,
    n_c = 30,
    n_t = 30,
    M = -0.1,
    threshold = 0.9,
    prior = "flat",
    prior_args = list(
      a_base_c = 1, b_base_c = 1,
      a_base_t = 1, b_base_t = 1
    ),
    n_draws = n_draws,
    show_progress = FALSE
  )

  # Call R version
  res_r <- bcts_power_betaBinom_conj(
    B = settings$B,
    p_c = settings$p_c,
    p_t = settings$p_t,
    n_c = settings$n_c,
    n_t = settings$n_t,
    M = settings$M,
    threshold = settings$threshold,
    prior = settings$prior,
    prior_args = settings$prior_args,
    n_draws = settings$n_draws,
    show_progress = settings$show_progress,
    seed = 42
  )

  # Call C++ version (returns logical vector of decisions)
  set.seed(42)
  res_cpp <- rct_power_beta_binom_cpp_vec(
    B = settings$B,
    p_c = settings$p_c,
    p_t = settings$p_t,
    n_c = settings$n_c,
    n_t = settings$n_t,
    M = settings$M,
    threshold = settings$threshold,
    prior = settings$prior,
    prior_args = settings$prior_args,
    n_draws = settings$n_draws,
    show_progress = settings$show_progress
  )

  # Compare estimated power (mean decision in cpp vs estimate in R)
  power_cpp <- mean(res_cpp)
  power_r   <- res_r$estimate

  expect_true(abs(power_cpp - power_r) < 0.01,
              info = paste("Power estimates differ too much:",
                           "R =", round(power_r, 4),
                           "C++ =", round(power_cpp, 4)))
})
