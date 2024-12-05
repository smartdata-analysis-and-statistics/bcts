
sim <- bcts(n_dose_sel = 60, n_ss_reest = 175, n_pln = 200, n_max = 250,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_ref = "Placebo", method = "mcmc", alpha = 0.01,
            th.fut = 0.15, th.eff = 0.9, th.prom = 0.5, ##
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/CIPO_01.rds")

sim <- bcts(n_dose_sel = 60, n_ss_reest = 150, n_pln = 200, n_max = 250,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_ref = "Placebo", method = "mcmc", alpha = 0.01,
            th.fut = 0.15, th.eff = 0.9, th.prom = 0.5, ##
            nsim = 10000, progress.bar = "text")
save(sim, file = "CIPO/CIPO_02.RData")

sim <- bcts(n_dose_sel = 60, n_ss_reest = 160, n_pln = 200, n_max = 250,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_ref = "Placebo", method = "mcmc", alpha = 0.01,
            th.fut = 0.15, th.eff = 0.9, th.prom = 0.5, ##
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/CIPO_02.rds")

sim <- bcts(n_dose_sel = 60, n_ss_reest = 170, n_pln = 200, n_max = 250,
            mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
            sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
            trt_ref = "Placebo", method = "mcmc", alpha = 0.01,
            th.fut = 0.15, th.eff = 0.9, th.prom = 0.5, ##
            nsim = 10000, progress.bar = "text")
saveRDS(sim, file = "CIPO/CIPO_03.rds")
