#' Plot of PPoS at Interim by Futility and Trial Outcome
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotInterimPPoS
plotInterimPPoS <- function(...) {
  UseMethod("plotInterimPPoS")
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

  # Calculate power without futility adjustment
  power_no_fut <- power(x, adjust_for_futility = FALSE, conf.level = conf.level)

  # Calculate power with futility adjustment
  power_with_fut <- power(x, adjust_for_futility = TRUE, conf.level = conf.level)

  fut.trig <- mc_error_proportion(x = sum(x$simresults$fut.trig),
                                  n = nrow(x$simresults),
                                  level = conf.level)

  inc.ss <- mc_error_proportion(x = sum(x$simresults$inc.ss),
                                n = nrow(x$simresults),
                                level = conf.level)

  # Initialize an empty data frame
  out <- data.frame(statistic = character(),
                    adjustment = character(),
                    est = numeric(),
                    cil = numeric(),
                    ciu = numeric())

  # Add power without futility adjustment
  out <- out %>% add_row(data.frame(statistic = "Power",
                                    adjustment = "Without Futility Adjustment",
                                    est = power_no_fut$est,
                                    cil = power_no_fut$lower,
                                    ciu = power_no_fut$upper))

  # Add power with futility adjustment
  out <- out %>% add_row(data.frame(statistic = "Power",
                                    adjustment = "With Futility Adjustment",
                                    est = power_with_fut$est,
                                    cil = power_with_fut$lower,
                                    ciu = power_with_fut$upper))

  # Add other probabilities
  out <- out %>% add_row(data.frame(statistic = "Pr(fut.trig)",
                                    adjustment = "N/A",
                                    est = fut.trig$est,
                                    cil = fut.trig$lower,
                                    ciu = fut.trig$upper))
  out <- out %>% add_row(data.frame(statistic = "Pr(inc.ss)",
                                    adjustment = "N/A",
                                    est = inc.ss$est,
                                    cil = inc.ss$lower,
                                    ciu = inc.ss$upper))

  # Convert to percentages
  out$est <- out$est * 100
  out$cil <- out$cil * 100
  out$ciu <- out$ciu * 100

  # Set the factor levels for 'adjustment' to control the order of the plot
  out$adjustment <- factor(out$adjustment,
                           levels = c("Without Futility Adjustment",
                                      "With Futility Adjustment",
                                      "N/A"))

  # Create the plot
  ggplot(out, aes(x = .data$statistic, y = .data$est, color = .data$adjustment)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = .data$cil, ymax = .data$ciu),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    facet_wrap(~statistic, scales = "free") +
    labs(x = "Statistic", y = "Estimate (95% CI)", color = "Adjustment") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

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

    color_values <- c("Dose 1" = "#3f88ab",
                      "Dose 2" = "#8e1f20",
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

  # Update decision labels with proportions in brackets
  facet_labels <- ben %>%
    dplyr::group_by(.data$decision) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prop = .data$n / sum(.data$n),
                  label = paste0(.data$decision, " (", round(.data$prop * 100, 1), "%)"))

  # Convert the labels to a named vector for facet labeling
  facet_labels_vec <- setNames(facet_labels$label, facet_labels$decision)


  x_label <- paste("Estimated benefit", dose_1)
  y_label <- paste("Estimated benefit", dose_2)

  # Plot the bar chart with percentage labels
  ggplot(ben, aes(x = ben[,1], y = ben[,2])) +
    geom_point(aes(color = .data$decision), size = 1, alpha = 0.5, show.legend = FALSE) +
    stat_density_2d(aes(fill = after_stat(.data$level)), geom = "polygon", alpha = 0.4) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    scale_fill_viridis_c(option = "plasma") +  # Color scale for density
    xlab(x_label) +
    ylab(y_label) +
    theme_minimal() +
    theme(legend.position = "right") +  # Move legend to the right
    guides(fill = guide_colorbar(title = "Density Level"), color = "none")  +# No color legend for decision, only for density level
    scale_color_manual(values = color_values) +
    facet_wrap(~.data$decision, labeller = as_labeller(facet_labels_vec))
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
    geom_point(size = 2, alpha = 0.5) +
    scale_color_gradientn(colors = c("#46125b", "#2b9689", "#f0e42e")) +  # Adjust colors as needed
    labs(x = "Estimated Benefit", y = "Predictive Power", color = "SS Increase") +
    theme_minimal() +
    theme(legend.position = "top",  # Move the legend to the top
          legend.direction = "horizontal")  # Arrange the legend items horizontally

}

#' Plot of PPoS at Interim by Futility and Trial Outcome
#'
#' @param x An object of class bcts
#' @param interim_index An integer specifying which interim analysis to use for the plot. If `NULL` (default), the latest interim analysis is used.
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @method plotInterimPPoS bcts
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr n
#' @importFrom rlang .data
plotInterimPPoS.bcts <- function(x, interim_index = NULL, ...) {

  if (is.null(interim_index)) {
    interim_index <- length(x$ss_inc)
  }

  inc.ss <- x$ss_inc[[interim_index]] > 0
  ppos <- x$PPos[[interim_index]] %>% pull("Interim")
  rejectH0 <- x$simresults$rejectH0.final
  fut.trig <- x$simresults$fut.trig

  ggdat <- data.frame(ppos = ppos, inc.ss = inc.ss,
                      fut.trig = fut.trig,
                      rejectH0 = rejectH0, decision = NA)
  ggdat$decision[which(ggdat$rejectH0)] <- "Trial Success"
  ggdat$decision[which(!ggdat$rejectH0)] <- "Trial Failure"
  ggdat$decision[which(ggdat$fut.trig)] <- "Futility"

  # Summarize to get the count and proportions for each group
  facet_labels <- ggdat %>%
    dplyr::group_by(fut.trig) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prop = .data$n / sum(.data$n),  # Calculate proportion
                  label = dplyr::case_when(
                    fut.trig == TRUE ~ paste0("Futility (", round(.data$prop * 100, 1), "%)"),
                    fut.trig == FALSE ~ paste0("Continue/Expand (", round(.data$prop * 100, 1), "%)")
                  ))

  # Convert the labels to a named vector for the facet labeller
  facet_labels_vec <- setNames(facet_labels$label, facet_labels$fut.trig)

  ggplot(ggdat, aes(x = .data$ppos, fill = .data$rejectH0)) +
    geom_histogram(binwidth = 0.02, boundary = 0, alpha = 0.85) +
    scale_fill_manual(values = c("TRUE" = "#66c2a5",  # Green for trial success
                                 "FALSE" = "#fc8d62"),  # Red for trial failure
                      labels = c("TRUE" = "Trial Success", "FALSE" = "Trial Failure")) +  # Updated labels
    xlab("Posterior Probability of Success (PPoS) for planned sample size") +
    labs(fill = "Trial Outcome") +
    ggtitle(paste("Distribution of PPoS at Interim", interim_index),
            subtitle = paste("Futility Analysis Across All Interim Looks")) +  # Added subtitle
    ylab("Number of simulations") +
    theme(legend.position = "top", legend.title = element_blank()) +
    facet_wrap(~ fut.trig, labeller = as_labeller(facet_labels_vec)) +
    theme_minimal()
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


#' Plot the Relationship Between Type-I Error and Power Across Decision Thresholds (Gamma)
#'
#' This function creates a scatter plot to visualize the relationship between Type-I error and power
#' across different decision thresholds (gamma values). The input is two lists of simulation results:
#' one for power and one for Type-I error. The resulting plot shows how power and Type-I error
#' vary with respect to gamma, and optionally distinguishes between different design criteria.
#'
#' @param sim_list_power A list of simulation results for power, where each element is an object containing power estimates.
#' @param sim_list_type1 A list of simulation results for Type-I error, where each element is an object containing Type-I error estimates.
#'
#' @details
#' The function expects the lists `sim_list_power` and `sim_list_type1` to be parallel, meaning each corresponding
#' element in the lists should refer to the same simulation setting (i.e., same gamma, number of looks, etc.).
#' It extracts the gamma values, power, Type-I error, and optionally other design-related parameters from
#' these lists and generates a scatter plot. Points in the plot are colored based on the gamma values and
#' shaped based on the prioritization of certain (e.g., lower) doses.
#'
#' @return A ggplot2 object representing the relationship between Type-I error and power across different decision thresholds (gamma).
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @export
plotPowerType1 <- function(sim_list_power, sim_list_type1) {

  # Create a data frame for plotting
  ggdat <- data.frame(gamma =  numeric(),
                      type1 = numeric(),
                      th.fut = numeric(),
                      th.eff = numeric(),
                      th.prom = numeric(),
                      power = numeric(),
                      no.looks = numeric(),
                      pld = logical())

  for (sim in sim_list_power) {
    result <- data.frame(gamma =  sim$gamma,
                         type1 = NA,
                         th.fut = sim$th.fut,
                         th.eff = sim$th.eff,
                         th.prom = sim$th.prom,
                         power = power(sim, adjust_for_futility = TRUE)$est,
                         no.looks = sim$no.looks,
                         pld = ifelse(is.null(sim$trt_rank), FALSE, TRUE))
    ggdat <- ggdat %>% add_row(result)
  }

  for (sim in sim_list_type1) {
    pld_i <- ifelse(is.null(sim$trt_rank), FALSE, TRUE)
    row <- which(ggdat$gamma == sim$gamma &
                   ggdat$th.fut == sim$th.fut &
                   ggdat$th.eff == sim$th.eff &
                   ggdat$th.prom == sim$th.prom &
                   ggdat$no.looks == sim$no.looks &
                   ggdat$pld == pld_i)
    ggdat$type1[row] <- power(sim, adjust_for_futility = FALSE)$est
  }
  ggdat$th.fut <- as.factor(ggdat$th.fut)

  # Plot Type-I error vs Power with gamma as color
  ggplot(ggdat, aes(x = .data$type1,
                    y = .data$power,
                    color = factor(.data$gamma),
                    shape = factor(.data$pld))) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      x = "Type-I Error",
      y = "Power",
      color = "Gamma",
      title = paste("Type-I Error vs Power")
    ) +
    #theme_minimal() +
    theme(legend.position = "right") +
    scale_color_viridis_d() +
    scale_shape_manual(values = c(16, 17),  # Different point shapes for pld = TRUE/FALSE
                       labels = c("Prioritize Most Effective Dose", "Prioritize Low Dose")) +
    facet_wrap(~ no.looks, labeller = as_labeller(c(`1` = "One Look",
                                                    `2` = "Two Looks")))
}

