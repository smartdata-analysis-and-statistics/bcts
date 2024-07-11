# Evaluate type-I error over a grid of gamma
nsim <- 50000
gammas <- seq(from = 0.98, to = 0.99, by = 0.001)

results <- data.frame(gamma = numeric(),
                      statistic = character(),
                      est = numeric(),
                      cil = numeric(),
                      ciu = numeric())

for (i in seq(gammas)) {
  sim <- bcts(n_int = 20*3,
              n_pln = 20 + 65*2,
              n_max = 20 + 80*2,
              mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5),
              sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1),
              trt_ref = "Placebo",
              prioritize_low_rank = TRUE, gamma = gammas[i], #   %
              method = "mcmc",
              th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
              nsim = 10000)
  results <- results %>% add_row(cbind(gamma = gammas[i], sim))
}

saveRDS(results, file = "CIPO/sim60.150.180.rds")

results2 <- data.frame(gamma = numeric(),
                      statistic = character(),
                      est = numeric(),
                      cil = numeric(),
                      ciu = numeric())

for (i in seq(gammas)) {
  sim <- bcts(n_int = 20*3,
              n_pln = 20 + 70*2,
              n_max = 20 + 80*2,
              mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
              sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
              trt_ref = "Placebo",
              prioritize_low_rank = TRUE, gamma = gammas[i], #   %
              method = "mcmc",
              th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
              nsim = 10000)
  results2 <- results2 %>% add_row(cbind(gamma = gammas[i], sim))
}

saveRDS(results2, file = "CIPO/sim20.70.80.rds")



sim98 <- bcts(n_int = 20*3,
                   n_pln = 20 + 65*2,
                   n_max = 20 + 80*2,
                   mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
                   sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
                   trt_ref = "Placebo",
                  prioritize_low_rank = TRUE, gamma = 0.98, #   %
                  method = "mcmc",
                  th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
                  nsim = 10000,
                  num_chains = 4,
                  n.iter = 10, #
                  n.adapt = 500,
                  perc_burnin = 0.2,
                  progress.bar = "text")
saveRDS(sim98, file = "CIPO/gamma98N20.65.80.rds")


sim98b <- bcts(n_int = 20*3,
              n_pln = 20 + 65*2,
              n_max = 20 + 80*2,
              mu = c("Placebo" = 0, "Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
              sigma = c("Placebo" = 1, "Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
              trt_ref = "Placebo",
              prioritize_low_rank = TRUE, gamma = 0.98, #   %
              method = "bayes",
              th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
              nsim = 10000,
              num_chains = 4,
              n.iter = 10, #
              n.adapt = 500,
              perc_burnin = 0.2,
              progress.bar = "text")
saveRDS(sim98, file = "CIPO/gamma98N20.65.80.rds")



