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
      p(HTML('For documentation and source code, visit the <a href="https://github.com/smartdata-analysis-and-statistics/bcts" target="_blank">GitHub repository</a>.')),
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

    ver <- as.character(utils::packageVersion("bcts"))
    date <- utils::packageDescription("bcts")$Date %||% Sys.Date()
    desc <- utils::packageDescription("bcts")

    author_raw <- desc$Author %||% "Unknown"

    # Clean up ORCID formatting (optional)
    author_clean <- gsub("\\s*\\(ORCID:[^\\)]*\\)", "", author_raw)  # remove ORCID line
    author_clean <- gsub("\\s+", " ", author_clean)  # remove excess whitespace

    output$version_text <- renderText({
      paste0(
        "bcts version: ", ver, " (", date, ")",
        "\nAuthor: ", author_clean,
        "\nR version: ", R.version.string,
        "\nPlatform: ", R.version$platform
      )
    })
  })
}
