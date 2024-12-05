#' Extract Treatment Names
#'
#' @param mu Named vector of expected effect sizes for each treatment group (including the control group).
#' @param sigma Named vector of expected standard deviations for each treatment group.
#' @return A named list with unique treatment names.
extract_treatment_names <- function(mu, sigma) {
  # Combine the names from both vectors
  combined_names <- union(names(mu), names(sigma))

  # If names are missing, assign default names
  if (is.null(combined_names) || any(combined_names == "")) {
    num_treatments <- max(length(mu), length(sigma))
    default_names <- paste("Treatment", seq_len(num_treatments))

    # Assign default names where missing
    if (is.null(names(mu)) || any(names(mu) == "")) {
      names(mu) <- default_names[seq_along(mu)]
    }
    if (is.null(names(sigma)) || any(names(sigma) == "")) {
      names(sigma) <- default_names[seq_along(sigma)]
    }

    # Recompute combined names
    combined_names <- union(names(mu), names(sigma))
  }

  return(combined_names)
}
