# Evaluate type-I error over a grid of gamma
require(binom)

nsim <- 1000
gammas <- seq(from = 0.98, to = 0.99, by = 0.001)

results <- data.frame(gamma = numeric(),
                      statistic = character(),
                      est = numeric(),
                      cil = numeric(),
                      ciu = numeric())

for (i in seq(gammas)) {
  sim <- eval_bcts(n_per_arm_int = 20,
                   n_per_arm_pln = 65,
                   n_per_arm_max = 80,
                   block.sizes = 1,
                   mu_c = 0, sd_c = 1,
                   mu_t = c("Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
                   sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
                   prioritize_low_trteff = TRUE,
                   gamma = gammas[i], #   %
                   th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
                   nsim = nsim,
                   num_chains = 4,
                   n.iter = 10, #
                   n.adapt = 500,
                   perc_burnin = 0.2,
                   progress.bar = "text", #
                   quiet = FALSE)
  results <- results %>% add_row(cbind(gamma = gammas[i], sim))
}


# Plot type-I error
ggplot(results %>% filter(statistic %in% c("Power", "Type-I error")), aes(x=gamma, y = est, ymin = cil, ymax = ciu, color = statistic)) +
  geom_point() +
  xlab("Gamma") +
  geom_errorbar() +
  geom_hline(yintercept = 0.03, lty = 2, col = "red") +
  geom_hline(yintercept = 0.025, lty = 2, col = "blue") +
  geom_hline(yintercept = 0.8, lty = 2, col = "grey") +
  ggtitle("N = 20; 65; 80")
