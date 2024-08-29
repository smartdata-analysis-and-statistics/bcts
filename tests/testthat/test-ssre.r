require(rjags)

test_that("Test sample size re-estimation at interim", {
  nsim <- 10000
  alpha <- 0.05

  mu <- c("Placebo" = 0, "Drug A" = 0.5)
  sigma <- c("Placebo" = 1, "Drug A" = 1)
  trt_ref <- "Placebo"
  n_pln <- 20 + 65*2
  n_max <- 20 + 80*2
  th.fut <- 0.2
  th.eff <- 0.9
  th.prom <- 0.5

  # Extract the treatment names
  trt_names <- names(mu)
  names(mu) <- names(sigma) <- trt_names

  # Active treatment names
  trt_active <- setdiff(trt_names, trt_ref)

  resultf  <- data.frame(fut.trig = rep(NA, nsim), inc.ss = NA,
                        n.final = NA, ppos_pln = NA, ppos_fin = NA, benefit = NA)

  for (i in 1:nsim) {

    # Simulate trial data at interim stage
    dat_int <- sim_rct_normal(n = 20*2,
                              mean = mu,
                              sd = sigma,
                              trtnames = trt_names)

    resultf[i,] <- ssre(dat_int = dat_int,
                              n_pln = n_pln,
                    n_max = n_max,
                              trt_ref = trt_ref, trt_active = trt_active,
                              gamma = 1 - alpha/2,
                              th.fut = th.fut, th.eff = th.eff,
                              method = "mcmc")$result
  }

  expect_false(any(resultf %>% filter(fut.trig) %>% pull(ppos_fin) > th.fut),
              info = "PPoS should be below th.fut when futility is triggered")




})
