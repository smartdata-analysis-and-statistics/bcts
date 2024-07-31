################################################################################
# Default scenario
################################################################################
n_int <- 20*3
n_pln <- 20 + 65*2
n_max <- 20 + 80*2
gamma <- 0.98
th.fut <- 0.2
th.eff <- 0.9
th.prom <- 0.5

### Frequentist
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = gamma, method = "mcmc",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom, ##
            nsim = 100000, progress.bar = "text")
saveRDS(sim, file = "CIPO/sim_type1_G98.N20.65.80.rds")

sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = gamma, method = "mcmc",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom, ##
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/sim_power_G98.N20.65.80.rds")

### Bayesian
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = gamma, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 100000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_type1_G98.N20.65.80.rds")

sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = gamma, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_power_G98.N20.65.80.rds")

## Gamma = 0.981
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = 0.981, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_type1_G981.N20.65.80.rds")

## Gamma = 0.982
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = 0.982, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_type1_G982.N20.65.80.rds")

## Gamma = 0.983
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = 0.983, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_type1_G983.N20.65.80.rds")
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = 0.983, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_power_G983.N20.65.80.rds")

## Gamma = 0.985
sim <- bcts(n_int = n_int, n_pln = n_pln, n_max = n_max,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_rank = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 2, "Placebo" = 3),
            trt_ref = "Placebo",
            prioritize_low_rank = TRUE, gamma = 0.985, method = "bayes",
            th.fut = th.fut, th.eff = th.eff, th.prom = th.prom,
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/simb_type1_G985.N20.65.80.rds")
