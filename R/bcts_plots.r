#' Plot Probability of trial success
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotPPoS
plotPPoS <- function(...) {
  UseMethod("plotPPoS")
}

#' Plot selected dose
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotSelDose
plotSelDose <- function(...) {
  UseMethod("plotSelDose")
}

#' Plot Interim Decision Outcomes
#'
#' This generic function creates a plot to visualize the interim decision outcomes from adaptive trial simulations.
#' The specific method used depends on the class of the input object.
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotInterimDecisions
plotInterimDecisions <- function(...) {
  UseMethod("plotInterimDecisions")
}

#' Plot Sample Size Re-Estimation Results
#'
#' This function plots the sample size re-estimation results from a simulation study. Each point on the plot
#' represents a simulation, with the x-axis showing the estimated benefit, the y-axis showing the predictive power,
#' and the color indicating the magnitude of the sample size increase.
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A ggplot2 object representing the plot of sample size re-estimation results.
#'
#' @details
#' The function extracts the estimated benefit and predictive power for the selected dose in each simulation and then
#' plots these against each other. The points are colored according to the magnitude of the sample size increase (`ss_inc`).
#' The color gradient legend is positioned at the top of the plot.
#'
#' @author Thomas Debray
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @export plotSampleSizeReEstimation
plotSampleSizeReEstimation <- function(...) {
  UseMethod("plotSampleSizeReEstimation")
}

#' Plot the final sample size
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotFinalSampleSize
plotFinalSampleSize  <- function(...) {
  UseMethod("plotFinalSampleSize")
}


#' Plot simulation results
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom rlang .data
plot.bcts <- function(x, ...) {
  conf.level <- 0.95

  power <- mc_error_proportion(x = sum(x$simresults$rejectH0.final),
                               n = nrow(x$simresults),
                               level = conf.level)

  fut.trig <- mc_error_proportion(x = sum(x$simresults$fut.trig),
                                  n = nrow(x$simresults),
                                  level = conf.level)

  inc.ss <- mc_error_proportion(x = sum(x$simresults$inc.ss),
                                n = nrow(x$simresults),
                                level = conf.level)

  out <- data.frame(statistic = character(),
                    est = numeric(),
                    cil = numeric(),
                    ciu = numeric())

  out <- out %>% add_row(data.frame(statistic = "Power",
                                    est = power$est,
                                    cil = power$lower,
                                    ciu = power$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(fut.trig)",
                                    est = fut.trig$est,
                                    cil = fut.trig$lower,
                                    ciu = fut.trig$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(inc.ss)",
                                    est = inc.ss$est,
                                    cil = inc.ss$lower,
                                    ciu = inc.ss$upper))

  # Convert to percentages
  out$est <- out$est * 100
  out$cil <- out$cil * 100
  out$ciu <- out$ciu * 100

  ggplot(out, aes(x = .data$statistic, y = .data$est)) +
    #geom_bar(stat = "identity", alpha = 0.7) +
    geom_point() +
    geom_errorbar(aes(ymin = .data$cil, ymax = .data$ciu), width = 0.2) +
    facet_wrap(~statistic, scales = "free") +
    labs(x = "Statistic", y = "Estimate (95% CI)") +
    xlab("") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) #+ scale_fill_brewer()

}

#' Plot posterior distribution of the mean outcome
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
posterior.bcts <- function(x, ...) {
  ggplot(x$simresults, aes(x = .data$est.final)) +
    geom_density() +
    xlab("Mean outcome selected dose")
}


#' Plot selected dose
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotSelDose.bcts <- function(x, ...) {

  # Calculate the percentage for each bar
  nDose <- x$simresults %>%
    dplyr::count(.data$sel.dose) %>%
    mutate(percentage = .data$n / sum(.data$n) * 100)

  # Plot the bar chart with percentage labels
  ggplot(nDose, aes(x = .data$sel.dose, y = .data$n, fill = .data$sel.dose)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(.data$percentage, 1), "%")), vjust = -0.5) +
    labs(x = "Selected Dose", y = "Count") +
    scale_fill_brewer()  +
    theme(legend.position = "none") +
    ylim(0,  nrow(x$simresults))

}

#' Interim Decision Outcomes: Comparative Effect Sizes by Dose
#'
#' Plot Interim Decision Outcomes for bcts Objects
#'
#' This method visualizes the interim decision outcomes from adaptive trial simulations for objects of class `bcts`.
#' It plots the estimated effect sizes of two different doses, with points colored according to the specified decision scheme.
#'
#' @param x An object of class `bcts`. This object should contain the simulation results and interim decisions related to dose selection and futility.
#' @param color_scheme A character string specifying the coloring scheme to be used. Options are "dose_selection" and "trial_decision".
#'                     The "dose_selection" scheme colors points based on dose selection decisions ("Select dose 1", "Select dose 2", "Futility").
#'                     The "trial_decision" scheme colors points based on trial decisions ("Continue", "Expand", "Futility").
#' @param interim_index An integer specifying which interim analysis to plot.
#'                      If `NULL` (default), the first interim analysis is used.
#' @param ... Optional arguments.
#' @return A `ggplot2` object showing the scatter plot of estimated effect sizes for the two doses, colored according to the specified decision scheme.
#'
#' @method plotInterimDecisions bcts
#' @export
#'
#' @seealso \code{\link{plotInterimDecisions}} for the generic function.
#' @author Thomas Debray
#'
#' @import ggplot2
#' @importFrom rlang .data
plotInterimDecisions.bcts <- function(x, color_scheme = "dose_selection",
                                      interim_index = NULL, ...) {

  if (!"benefit" %in% names(x)) {
    stop("No interim results are available!")
  }
  if (is.null(interim_index)) {
    interim_index <- 1
  }

  dose_1 <- colnames(x$benefit[[interim_index]])[1]
  dose_2 <- colnames(x$benefit[[interim_index]])[2]

  ben <- x$benefit[[interim_index]]
  ben$decision <- NA

  if (color_scheme == "dose_selection") {
    ben$decision[which(x$simresults$sel.dose == dose_1)] <- paste("Select", dose_1)
    ben$decision[which(x$simresults$sel.dose == dose_2)] <- paste("Select", dose_2)
    ben$decision[which(x$simresults$fut.trig)] <- "Futility"

    color_values <- c("Select Velusetrag 15mg" = "#3f88ab",
                      "Select Velusetrag 30mg" = "#8e1f20",
                      "Futility" = "#59b559")

  } else if (color_scheme == "trial_decision") {
    ben$decision[which(!x$simresults$inc.ss)] <- "Continue"
    ben$decision[which(x$simresults$inc.ss)] <- "Expand"
    ben$decision[which(x$simresults$fut.trig)] <- "Futility"

    color_values <- c("Continue" = "#59b559",
                      "Expand" = "#cda026",
                      "Futility" = "#8e1f20")
  } else {
    stop("Invalid color_scheme. Choose 'dose_selection' or 'trial_decision'.")
  }

  x_label <- paste("Estimated benefit", dose_1)
  y_label <- paste("Estimated benefit", dose_2)

  # Plot the bar chart with percentage labels
  ggplot(ben, aes(x = ben[,1], y = ben[,2], color = .data$decision)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    xlab(x_label) +
    ylab(y_label) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title = NULL)) +
    scale_color_manual(values = color_values)
}

#' Plot Sample Size Re-Estimation Results
#'
#' This function plots the sample size re-estimation results from a simulation study. Each point on the plot
#' represents a simulation, with the x-axis showing the estimated benefit, the y-axis showing the predictive power,
#' and the color indicating the magnitude of the sample size increase.
#'
#' @param x A list object containing simulation results, including benefit estimates, predictive power, and sample size increases.
#'          The list must include elements named `benefit`, `PPos`, and `ss_inc`, each of which should be a list containing matrices or data frames.
#' @param interim_index An integer specifying which interim analysis to use for the plot. If `NULL` (default), the latest interim analysis is used.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A ggplot2 object representing the plot of sample size re-estimation results.
#'
#' @details
#' The function extracts the estimated benefit and predictive power for the selected dose in each simulation and then
#' plots these against each other. The points are colored according to the magnitude of the sample size increase (`ss_inc`).
#' The color gradient legend is positioned at the top of the plot.
#'
#' @author Thomas Debray
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @method plotSampleSizeReEstimation bcts
#' @export
plotSampleSizeReEstimation.bcts <- function(x, interim_index = NULL, ...) {
  if (!"benefit" %in% names(x)) {
    stop("No interim results are available!")
  }
  if (is.null(interim_index)) {
    interim_index <- length(x$benefit)
  }

  ben <- x$benefit[[interim_index]]
  ppos <- x$PPos[[interim_index]]

  # Initialize an empty vector to store the benefit values
  selected_benefits <- selected_ppos <- numeric(length(x$simresults$sel.dose))

  # Loop through each simulation and extract the corresponding benefit
  for (i in seq_along(x$simresults$sel.dose)) {
    selected_dose <- x$simresults$sel.dose[i]  # Get the selected dose for this simulation
    selected_benefits[i] <- ben[i, selected_dose]  # Extract the corresponding benefit
    selected_ppos[i] <- ppos[i, "Interim"]
  }

  ggdat <- data.frame(ben = selected_benefits, ppos = selected_ppos,
                      ss_inc = x$ss_inc[[interim_index]])
  ggplot(ggdat, aes(x = .data$ben, y = .data$ppos, color = .data$ss_inc)) +
    geom_point() +
    scale_color_gradient(low = "#46125b", high = "#f0e42e") +  # Adjust colors as needed
    labs(x = "Estimated Benefit", y = "Predictive Power", color = "SS Increase") +
    theme_minimal() +
    theme(legend.position = "top",  # Move the legend to the top
          legend.direction = "horizontal")  # Arrange the legend items horizontally

}

#' Plot Probability of trial success
#'
#' @param x An object of class bcts
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @method plotPPoS bcts
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotPPoS.bcts <- function(x, ...) {

  pr.ss.increase <- round(mean(x$simresults$inc.ss)*100,0)
  inc.ss.labs <- c(paste0("Planned sample size increased (", pr.ss.increase, "%)"),
                   paste0("Planned sample size retained (", 100 - pr.ss.increase, "%)"))
  names(inc.ss.labs) <- c("TRUE", "FALSE")

  ggplot(x$simresults, aes(x = .data$sel.dose.ppos, fill = .data$rejectH0.final)) +
    geom_histogram(binwidth = 0.025, boundary = 0) +
    scale_fill_manual(values = c("FALSE" = "#B8B8B8", "TRUE" = "#1A80BB"),
                      labels = c("FALSE" = "Trial failure", "TRUE" = "Trial success")) +
    xlab("PPoS at final sample size") +
    labs(fill = "Trial success") +
    theme(legend.position = "bottom") +
    ylab("Number of simulations") +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_grid(~ inc.ss, labeller = labeller(inc.ss = inc.ss.labs), scales = "free_x")
}

#' Plot final sample size
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotFinalSampleSize.bcts <- function(x, ...) {

  # Join n.final with selected treatment
  #dat <- x$n.final %>% merge(x$simresults %>% select("sim", "sel.dose"))

  #ggplot(dat, aes(x=n, group = Treatment, fill = Treatment)) +
  #  geom_histogram(col = "black") +
  #  facet_wrap(~sel.dose) + scale_fill_brewer()

  # Calculate the mean final sample size
  mean_n_final <- mean(x$simresults$n.final)

  ggplot(x$simresults, aes(x = .data$n.final)) +
    geom_histogram() +
    xlab("Final sample size") +
    ylab("Count") +
    annotate("text", x = Inf, y = Inf, label = paste("Mean final sample size:", round(mean_n_final, 0)),
             hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +
    ylim(0,  nrow(x$simresults))
}
