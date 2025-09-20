#' Run the BCTS Shiny App
#'
#' Launches the interactive Shiny app included in the \code{bcts} package, which supports
#' Bayesian clinical trial simulations for both single-arm and randomized trials. This app provides
#' user-friendly controls to explore design parameters, power, and Type-I error.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}, such as host or port.
#' @param launch.browser Logical. Should the app be launched in the system browser? Defaults to \code{TRUE}.
#'
#' @details This function locates the Shiny app bundled in the \code{inst/app} directory of the
#' \code{bcts} package and runs it using \code{shiny::runApp}. Ensure that the package was installed
#' with the app directory included.
#'
#' @return None. This function is called for its side effect of launching the app.
#'
#' @examples
#' if (interactive()) {
#'   run_bcts_app()
#' }
#'
#' @export
run_bcts_app <- function(..., launch.browser = TRUE) {
  app_dir <- system.file("app", package = "bcts")
  if (app_dir == "") {
    stop("App not found. Try reinstalling the bcts package.", call. = FALSE)
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
