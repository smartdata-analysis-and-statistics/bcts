test_that("Analytic and simulated Type-I error are consistent for single-arm Beta-Binomial design", {
  # Parameters
  n_t       <- 35
  M         <- 0.6
  threshold <- 0.90
  a_base    <- 1
  b_base    <- 1
  p_null    <- M     # Evaluate at the boundary

  # Exact Type-I error (no Monte Carlo)
  exact_res <- .Call(
    `_bcts_sat_betabinom_type1_exact`,
    n_t, M, threshold,
    a_base, b_base, p_null
  )
  exact_type1 <- exact_res$estimate

  # Simulated Type-I error
  set.seed(123)
  sim_res <- .Call(
    `_bcts_sat_betabinom_type1`,
    50000,      # B
    n_t, M, threshold,
    a_base, b_base,
    FALSE       # show_progress
  )
  sim_type1 <- sim_res$estimate
  sim_se    <- sim_res$mc_se

  # Allow difference within 3 * Monte Carlo SE
  diff <- abs(sim_type1 - exact_type1)
  expect_true(
    diff < 3 * sim_se,
    sprintf(
      "Simulated type-I %.4f vs exact %.4f exceeds 3*mc_se = %.4f",
      sim_type1, exact_type1, 3 * sim_se
    )
  )
})
