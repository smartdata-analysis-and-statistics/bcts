#' Internal utilities for conjugate backend
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

.match_val <- function(x, y, tol = 1e-8) {
  if (is.null(y)) return(TRUE)
  if (is.numeric(x) && is.numeric(y)) return(abs(x - y) <= tol)
  identical(x, y)
}

# Build (alpha,beta) prior for one arm under flat or power prior
# prior_kind: "flat" or "power"
# args list: list(a_base = 1, b_base = 1, y0, n0, a0)
.conj_prior_shapes <- function(prior_kind = c("flat","power"), args = list()) {
  prior_kind <- match.arg(prior_kind)
  a_base <- args$a_base %||% 1
  b_base <- args$b_base %||% 1
  if (prior_kind == "flat") {
    c(alpha = a_base, beta = b_base)
  } else {
    y0 <- args$y0 %||% stop("prior_args$y0 required for power prior")
    n0 <- args$n0 %||% stop("prior_args$n0 required for power prior")
    a0 <- args$a0 %||% 0.5
    c(alpha = a_base + a0 * y0,
      beta  = b_base + a0 * (n0 - y0))
  }
}

# Given prior shapes and binomial data, return posterior shapes
# prior is c(alpha=?, beta=?)
# returns c(alpha=?, beta=?)
.conj_posterior_shapes <- function(prior, y, n) {
  c(alpha = prior["alpha"] + y,
    beta  = prior["beta"]  + (n - y))
}
