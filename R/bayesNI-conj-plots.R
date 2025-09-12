library(ggplot2)
library(dplyr)

plot_prior_likelihood_weight <- function(res) {
  stopifnot(inherits(res, "bayesNI"))
  s  <- res$settings
  pa <- s$prior_args
  
  if (!identical(s$prior, "power"))
    stop("This visualization is only meaningful for power prior runs.")
  
  # Effective prior sample size
  ESS_prior <- pa$a0 * pa$n_0
  ESS_data  <- s$n_c
  
  df <- tibble(
    source = c("Power prior (ESS)", "Current data (ESS)"),
    ESS    = c(ESS_prior, ESS_data)
  )
  
  ggplot(df, aes(x = source, y = ESS, fill = source)) +
    geom_col(width = 0.6, color = "black") +
    geom_text(aes(label = round(ESS,1)), vjust = -0.5) +
    labs(
      title = "Relative weight of prior vs current control likelihood",
      subtitle = sprintf("Power prior: a0 = %.2f, n0 = %d (ESS = %.1f)\nCurrent control n = %d",
                         pa$a0, pa$n_0, ESS_prior, s$n_c),
      y = "Effective sample size",
      x = NULL
    ) +
    scale_fill_manual(values = c("Power prior (ESS)" = "steelblue", 
                                 "Current data (ESS)" = "orange")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}