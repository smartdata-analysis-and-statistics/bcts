
# We first evaluate the type-I error
sim1 <- eval_power(n_int = 60, n_pln = 150, n_max = 180, block_size = 6,
                   mu_c = 0, sd_c = 1,
                   mu_t = c("Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
                   sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
                   prioritize_low_trteff = TRUE,
                   gamma = 0.98, #   %
                   th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
                   nsim = 1000, # no. of simulatoins at trial start
                   nsim.ppos = 0, # no. of simulations to assess PPOS
                   nsim_int = 100, # no.of simulations at interim
                   num_chains = 4,
                   n.iter = 5000, #
                   n.adapt = 500,
                   perc_burnin = 0.2,
                   progress.bar = "text", #
                   quiet = FALSE)
saveRDS(sim1, file = "CIPO/sim1.rds")


#Draper D, Hodges JS, Mallows CL, Pregibon D. Exchangeability and data analysis. J Royal Stat Soc Series A, (Stat Soc). 1993;156(1):9-28.
