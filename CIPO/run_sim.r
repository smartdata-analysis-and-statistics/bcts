
# We first evaluate the type-I error
sim1 <- bcts(n_per_arm_int = 20,
             n_per_arm_pln = 65,
             n_per_arm_max = 80,
             block.sizes = 1,
                   mu_c = 0, sd_c = 1,
                   mu_t = c("Velusetrag 15mg" = 0, "Velusetrag 30mg" = 0), #
                   sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
                   prioritize_low_trteff = TRUE,
                   gamma = 0.98, #   %
                   th.fut = 0.2, th.eff = 0.9, th.prom = 0.5, ##
                   nsim = 5000, # no. of simulatoins at trial start
                   num_chains = 4,
                   n.iter = 1000, #
                   n.adapt = 500,
                   perc_burnin = 0.2,
                   progress.bar = "text", #
                   quiet = FALSE)
saveRDS(sim1, file = "CIPO/sim1.rds")

### Power for 15mg at interim stage
mean(sim1$results_interim %>% filter(trt==2) %>% pull(sigb))
mean(sim1$results_interim %>% filter(trt==3) %>% pull(sigb))

### Power for 15mg at max sample size



## TODO: asssess type 1 at 30mg , assess type 1 at 15 mg, and at selected dose


#Draper D, Hodges JS, Mallows CL, Pregibon D. Exchangeability and data analysis. J Royal Stat Soc Series A, (Stat Soc). 1993;156(1):9-28.
