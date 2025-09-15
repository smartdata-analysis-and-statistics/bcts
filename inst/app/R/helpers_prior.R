prior_args_from_inputs <- function(prior, a0, y0, n0, abase, bbase) {
  if (identical(prior, "flat")) return(list())
  list(a0 = a0/100, y_0 = y0, n_0 = n0, a_base = abase, b_base = bbase)
}
