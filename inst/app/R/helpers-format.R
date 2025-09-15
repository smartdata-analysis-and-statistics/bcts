fmt_pct  <- function(x, d=1) sprintf(paste0("%.", d, "f%%"), 100*x)
fmt_ci   <- function(lo, hi, d=1) sprintf(paste0("[%.", d, "f, %.", d, "f]%%"), 100*lo, 100*hi)
fmt_int  <- function(x) formatC(as.integer(x), big.mark=",", format="d")
fmt_prob <- function(x, d=3) sprintf(paste0("%.", d, "f"), x)
