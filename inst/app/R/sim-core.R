# R/sim-core.R
run_all_oc <- function(pt, nt, pc, nc, M, mode, gamma, alpha, calibrate_on,
                       prior, prior_args, B, ndraws, seed,
                       progress_fun = NULL, verbose = FALSE) {

  # Optional progress callback wrapper
  pf <- progress_fun
  if (!is.null(pf) && !is.function(pf)) pf <- NULL

  if (mode == "alpha") {
    cal <- bcts_calibrate_betaBinom_conj(
      alpha = alpha, p_c = pc, M = M, n_c = nc, n_t = nt,
      prior = prior, prior_args = prior_args,
      B_cal = B, n_draws = ndraws, seed = seed,
      show_progress = FALSE, verbose = verbose,
      progress_fun = pf,
      calibrate_on = calibrate_on
    )
    gamma_used <- cal$gamma
    cal_obj    <- cal
  } else {
    gamma_used <- gamma
    cal_obj    <- NULL
  }

  t1 <- bcts_type1_betaBinom_conj(
    B = B, p_c = pc, M = M, n_c = nc, n_t = nt,
    threshold = gamma_used, prior = prior, prior_args = prior_args,
    n_draws = ndraws, show_progress = FALSE
  )

  pw <- bcts_power_betaBinom_conj(
    B = B, p_c = pc, p_t = pt, n_c = nc, n_t = nt, M = M,
    threshold = gamma_used, prior = prior, prior_args = prior_args,
    n_draws = ndraws, seed = seed, show_progress = FALSE
  )

  list(
    decision_mode = mode,
    alpha_target  = if (mode == "alpha") alpha else NA_real_,
    gamma_used    = gamma_used,
    cal           = cal_obj,
    t1            = t1,
    pw            = pw
  )
}
