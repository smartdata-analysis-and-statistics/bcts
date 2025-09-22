test_that("Analytic and simulated power are consistent for single-arm Beta-Binomial design", {
  # Set parameters
  p_t       <- 0.75
  n_t       <- 35
  M         <- 0.6
  threshold <- 0.90
  prior     <- "flat"
  a_base    <- 1
  b_base    <- 1

  # Analytic power
  exact_res <- .Call(
    `_bcts_singlearm_beta_power_exact`,
    p_t,        # true response rate
    n_t,        # sample size
    M,          # decision threshold
    threshold,  # posterior probability threshold
    prior,      # prior type
    a_base,     # beta prior alpha
    b_base      # beta prior beta
  )
  exact_power <- exact_res$estimate

  # Simulated power (use large B for stability)
  set.seed(123)
  sim_res <- .Call(
    `_bcts_singlearm_beta_power`,
    50000,          # B
    p_t,            # true response rate
    n_t,            # sample size
    M,              # decision threshold
    threshold,      # posterior probability threshold
    prior,          # prior type
    a_base,         # beta prior alpha
    b_base,         # beta prior beta
    FALSE           # show_progress
  )

  sim_power <- sim_res$estimate
  sim_se <- sim_res$mc_se

  # Check that simulated estimate is within 3 * mc_se of the exact value
  diff <- abs(sim_power - exact_power)
  expect_true(
    diff < 3 * sim_se,
    sprintf(
      "Simulated power %.4f vs exact %.4f exceeds 3*mc_se = %.4f",
      sim_power, exact_power, 3 * sim_se
    )
  )
})
