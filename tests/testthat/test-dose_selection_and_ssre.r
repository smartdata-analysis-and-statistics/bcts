require(rjags)

test_that("Test adaptive dose selection and sample size re-estimation", {
  nsim <- 10000
  alpha <- 0.05

  mu <- c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0)
  sigma <- c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1)
  trt_rank <- c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3)
  trt_ref <- "Placebo"
  n_pln <- 20 + 65*2
  n_max <- 20 + 80*2
  gamma <- 0.98
  prioritize_low_rank <- TRUE
  th.fut <- 0.2
  th.eff <- 0.9
  th.prom <- 0.5

  # Extract the treatment names
  trt_names <- names(mu)
  names(mu) <- names(sigma) <- trt_names

  # Active treatment names
  trt_active <- setdiff(trt_names, trt_ref)

  sigf <- sigb <- rep(NA, nsim)

  inc.ss <- fut.trig  <- data.frame(freq = logical(), bayes = logical())
  n.final <- data.frame(freq = numeric(), bayes = numeric())

  for (i in 1:nsim) {

    # Simulate trial data at interim stage
    dat_int <- sim_rct_normal(n = 20*3,
                              mean = mu,
                              sd = sigma,
                              trtnames = trt_names)

    resultf <- dose_selection_and_ssre(dat_int = dat_int,
                                   n_pln = n_pln, n_max = n_max,
                                   trt_ref = trt_ref, trt_active = trt_active,
                                   trt_rank = trt_rank,
                                   prioritize_low_rank = prioritize_low_rank, gamma = gamma,
                                   th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
                                   method = "mcmc")

    resultb <- dose_selection_and_ssre(dat_int = dat_int,
                                   n_pln = n_pln, n_max = n_max,
                                   trt_ref = trt_ref, trt_active = trt_active,
                                   trt_rank = trt_rank,
                                   prioritize_low_rank = prioritize_low_rank, gamma = gamma,
                                   th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
                                   method = "bayes",
                                   num_chains = 4, n.iter = 5000, n.adapt = 500, perc_burnin = 0.2)

    inc.ss <- inc.ss %>% add_row(data.frame(freq = resultf$result$inc.ss, bayes = resultb$result$inc.ss))
    fut.trig <- fut.trig %>% add_row(data.frame(freq = resultf$result$fut.trig, bayes = resultb$result$fut.trig))
    n.final <- n.final %>% add_row(data.frame(freq = resultf$result$n.final, bayes = resultb$result$n.final))
  }

  test.fut.trig <- prop.test(x = colSums(fut.trig), n = rep(nsim, 2), alternative = "two.sided", correct = TRUE)
  expect_true(test.fut.trig$p.value > alpha,
              info = paste0("Futility triggering is not equally likely across the estimation methods (",
                            paste(colMeans(fut.trig), collapse = "; "), ")"))

  test.inc.ss <- prop.test(x = colSums(inc.ss), n = rep(nsim, 2), alternative = "two.sided", correct = TRUE)
  expect_true(test.inc.ss$p.value > alpha,
              info = paste0("Sample size increase is not equally likely across the estimation methods (",
                            paste(colMeans(inc.ss), collapse = "; "), ")"))

  test.n.final <- t.test(n.final$freq, n.final$bayes, var.equal = FALSE)
  expect_true(test.n.final$p.value > alpha,
              info = paste0("Final sample size is not equal across the estimation methods (",
                            paste(colMeans(n.final), collapse = "; "), ")"))

})
