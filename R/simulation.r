#' Simulate a randomized clinical trial
#'
#' @param n Sample size
#' @param mean vector with the means of the treatment groups
#' @param sd vector with the standard deviations of the treatment groups
#' @param trtnames vector with treatment names
#'
#' @return A dataframe with the simulated trial data
#'
#' @export
#' @importFrom stats rnorm
#' @importFrom dplyr mutate rename slice_head
#' @importFrom rlang .data
#' @importFrom randomizeR crPar genSeq getRandList
#'
sim_rct_normal <- function(n,
                           mean,
                           sd,
                           trtnames) {

  # Generate randomization parameters
  params <- crPar(N = n,
                  K = length(trtnames),
                  ratio = rep(1, length(trtnames)),
                  groups = trtnames)

  # Get randomization sequence
  rs <- genSeq(params)

  dat <- data.frame(id = seq(n),
                    Treatment = getRandList(rs)[1,],
                    Y = NA
                    )


  # Create a data frame to store the simulated trial data
  for (trt in seq(trtnames)) {
    ntrt <- sum(dat$Treatment == trtnames[trt])
    dat$Y[which(dat$Treatment == trtnames[trt])] <- rnorm(n = ntrt, mean = mean[trt], sd[trt])
  }
  dat <- dat %>% mutate(Treatment = factor(.data$Treatment))

  # Select n top records
  return(dat)
}

