require(rjags)
require(testthat)

test_that("Test RCT with one interim analysis", {

  sim <- bcts_one_interim(n_int = 60,
                          n_pln = 150,
                          n_max = 180,
                          mu =  c("Placebo" = 0, "Drug A" = 0.4, "Drug B" = 0.5),
                          sigma = c("Placebo" = 1, "Drug A" = 1, "Drug B" = 1),
                          trt_ref = "Placebo",
                          gamma = 0.975,
                          th.fut = 0.2, th.eff = 0.9, th.prom = 0.5,
                          method = "mcmc", nsim = 1000)

  ppos_int <- sim$PPos$Int1$Interim

  ind_ss_inc <- which(ppos_int >= 0.5 & ppos_int < 0.9)4
  ind_no_ss_inc <- which(ppos_int < 0.5 | ppos_int >= 0.9)

  expect_false(any(sim$ss_inc$Int1[ind_ss_inc] == 0))
  expect_false(any(sim$ss_inc$Int1[ind_no_ss_inc] > 0))

})
