#' Plot of PPoS at Interim by Futility and Trial Outcome

#' @param x An object of class bcts_results
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotInterimPPoS.bcts_results <- function(x, ...) {
  plotInterimPPoS(x$sim_power)
}


#' Plot final sample size
#'
#' @param x An object of class bcts_results
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotFinalSampleSize.bcts_results <- function(x,...) {
  plotFinalSampleSize(x$sim_power)
}
