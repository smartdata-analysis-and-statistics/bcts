#' @title Info UI Module
#' @description UI for version and basic info
#' @export
mod_info_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    title = "Info",
    fluidPage(
      h4("About this app"),
      p("This application allows simulation and evaluation of Bayesian trial designs using conjugate Beta–Binomial models."),
      p("It supports both randomized trials and single-arm studies with optional historical borrowing."),
      br(),
      h5("Version info"),
      verbatimTextOutput(ns("version_text"))
    )
  )
}

#' @title Info Server Module
#' @description Displays package version info
#' @export
mod_info_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$version_text <- renderText({
      paste(
        "bcts version:", as.character(utils::packageVersion("bcts")),
        "\nR version:", R.version.string,
        "\nPlatform:", R.version$platform
      )
    })
  })
}
