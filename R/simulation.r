sim_rct_normal <- function(n,
                           mean,
                           sd,
                           trtnames,
                           block.sizes) {

  # Treatment allocation
  dat <- blockrand(n = n,
                    num.levels = length(trtnames), # three treatments
                    levels = trtnames, # arm names
                    block.sizes = block.sizes) # stratum abbrev

  dat <- dat %>% mutate(Y = NA)

  # Create a data frame to store the simulated trial data
  for (trt in seq(trtnames)) {
    ntrt <- sum(dat$treatment == trtnames[trt])
    dat$Y[which(dat$treatment == trtnames[trt])] <- rnorm(n = ntrt, mean = mean[trt], sd[trt])
  }
  dat <- dat %>% mutate(treatment = factor(treatment),
                        trt = factor(treatment,
                                     labels = seq(trtnames),
                                     levels = trtnames)) %>%
    rename(Treatment = treatment)

  # Select n top records
  return(dat %>% slice_head(n = n))
}

eval_bcts <- function(n_per_arm_int = 20,
                      n_per_arm_pln = 65,
                      n_per_arm_max = 80,
                      block.sizes = 1, # block size for randomization
                   mu_c = 0,
                   sd_c = 1,
                   mu_t = c("Velusetrag 15mg" = 0.4, "Velusetrag 30mg" = 0.5), #
                   sd_t = c("Velusetrag 15mg" = 1, "Velusetrag 30mg" = 1), #
                   prioritize_low_trteff = TRUE,
                   gamma = 0.98, #   %
                   method = "mcmc",
                   th.fut = 0.2, ##
                   th.eff = 0.9, ##
                   th.prom = 0.5, ##
                   nsim = 1000, # no. of simulatoins at trial start
                   num_chains = 4,
                   n.iter = 5000, #
                   n.adapt = 500,
                   perc_burnin = 0.2, # How many n.iter should be reserved for burnin?
                   slevel = 0.05, # Confidence level for summary estimates of interest
                   progress.bar = "text") {

  # Evaluate type-I error
  mu_t_type1 <- mu_t
  mu_t_type1[] <- 0
  sim_type1 <- bcts(n_per_arm_int = n_per_arm_int, n_per_arm_pln = n_per_arm_pln,
                    n_per_arm_max = n_per_arm_max, block.sizes = block.sizes,
                    mu_c = 0, sd_c = sd_c, mu_t = mu_t_type1, sd_t = sd_t, #
                    prioritize_low_trteff = prioritize_low_trteff, gamma = gamma, #   %
                    method = method,
                    th.fut = th.fut, th.eff = th.eff, th.prom = th.prom, ##
                    nsim = nsim, num_chains = num_chains, n.iter = n.iter, #
                    n.adapt = n.adapt, perc_burnin = perc_burnin, progress.bar = progress.bar)

  # Evaluate power
  sim_power <- bcts(n_per_arm_int = n_per_arm_int, n_per_arm_pln = n_per_arm_pln,
                    n_per_arm_max = n_per_arm_max, block.sizes = block.sizes,
                    mu_c = mu_c, sd_c = sd_c, mu_t = mu_t, sd_t = sd_t, #
                    prioritize_low_trteff = prioritize_low_trteff, gamma = gamma, #   %
                    method = method,
                    th.fut = th.fut, th.eff = th.eff, th.prom = th.prom, ##
                    nsim = nsim, num_chains = num_chains, n.iter = n.iter, #
                    n.adapt = n.adapt, perc_burnin = perc_burnin, progress.bar = progress.bar)

  ## Get quantiles of final sample size
  n_fin_qt <- stats::quantile(sim_power$simresults$n_per_arm_fin, c(slevel/2, 1 - (slevel/2)))

  n_arm_dropped <- length(mu_t) - 1

  out <- data.frame(statistic = character(), est = numeric(), cil = numeric(), ciu = numeric())

  out <- out %>% add_row(data.frame(statistic = "Type-I error", est = sim_type1$power$mean, cil = sim_type1$power$lower, ciu = sim_type1$power$upper))
  out <- out %>% add_row(data.frame(statistic = "Power", est = sim_power$power$mean, cil = sim_power$power$lower, ciu = sim_power$power$upper))
  out <- out %>% add_row(data.frame(statistic = "N interim", est = sim_power$n$interim, cil = sim_power$n$interim, ciu = sim_power$n$interim))
  out <- out %>% add_row(data.frame(statistic = "N final",
                                    est = sim_power$n$final,
                                    cil = n_fin_qt[1]*2+n_arm_dropped*n_per_arm_int,
                                    ciu = n_fin_qt[2]*2+n_arm_dropped*n_per_arm_int))
  out <- out %>% add_row(data.frame(statistic = "Pr(fut.trig)",
                                    est = sim_power$fut.trig$mean,
                                    cil = sim_power$fut.trig$lower,
                                    ciu = sim_power$fut.trig$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(inc.ss)",
                                    est = sim_power$inc.ss$mean,
                                    cil = sim_power$inc.ss$lower,
                                    ciu = sim_power$inc.ss$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(sel.low)",
                                    est = sim_power$sel.low$mean,
                                    cil = sim_power$sel.low$lower,
                                    ciu = sim_power$sel.low$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(sel.high)",
                                    est = sim_power$sel.high$mean,
                                    cil = sim_power$sel.high$lower,
                                    ciu = sim_power$sel.high$upper))

  return(out)
}

blockrand <- function(n, num.levels = 2, levels = LETTERS[seq(length = num.levels)],
          id.prefix, stratum, block.sizes = 1:4, block.prefix, uneq.beg = FALSE,
          uneq.mid = FALSE, uneq.min = 0, uneq.maxit = 10)
{
  treat <- vector(mode(levels))
  block.id <- numeric(0)
  block.size <- numeric(0)
  i <- 1
  while (length(treat) < n) {
    if ((uneq.beg && (i == 1)) || (uneq.mid && (length(treat) > n/2))) {
      n2 <- max(block.sizes) * length(levels) + 1
      ii <- uneq.maxit
      while (ii > 0) {
        tmp <- sample(levels, n2, replace = TRUE)
        if (diff(range(table(tmp))) >= uneq.min) {
          ii <- 0
        }
        ii <- ii - 1
      }
      treat <- c(treat, tmp)
      block.id <- c(block.id, rep(i, n2))
      block.size <- c(block.size, rep(n2, n2))
      i <- i + 1
      if (length(treat) > n/2) {
        uneq.mid <- FALSE
      }
    }
    block.n <- if (length(block.sizes) > 1) {
      sample(block.sizes, 1)
    }
    else {
      block.sizes
    }
    tmp <- rep(levels, block.n)
    n2 <- length(tmp)
    treat <- c(treat, sample(tmp))
    block.id <- c(block.id, rep(i, n2))
    block.size <- c(block.size, rep(n2, n2))
    i <- i + 1
  }
  n3 <- length(treat)
  if (missing(id.prefix)) {
    id <- seq(length.out = n3)
  }
  else if (is.numeric(id.prefix)) {
    id <- seq(id.prefix, length.out = n3)
  }
  else {
    w <- floor(logb(n3, 10)) + 1
    id <- paste(id.prefix, formatC(seq(length.out = n3),
                                   flag = "0", width = w), sep = "")
  }
  out <- data.frame(id = I(id))
  if (!missing(stratum)) {
    out$stratum <- rep(stratum, n3)
  }
  if (missing(block.prefix)) {
    out$block.id <- factor(block.id)
  }
  else if (is.numeric(block.prefix)) {
    out$block.id <- factor(block.id + block.prefix)
  }
  else {
    w <- floor(logb(max(block.id), 10)) + 1
    out$block.id <- factor(paste(block.prefix, formatC(block.id,
                                                       flag = "0", width = w), sep = ""))
  }
  out$block.size <- block.size
  out$treatment <- factor(treat)
  return(out)
}
