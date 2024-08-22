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

#' Plot Interim Decisions
#'
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot2} object.
#'
#' @details This is a generic function.
#'
#' @export plotSelDose
plotInterimDecisions <- function(...) {
  UseMethod("plotInterimDecisions")
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

#' Plot interim decisions
#'
#' @param x An object of class bcts
#' @param ... Optional arguments
#'
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plotInterimDecisions.bcts <- function(x, ...) {

  dose_1 <- colnames(x$interim$benefit)[1]
  dose_2 <- colnames(x$interim$benefit)[2]

  ben <- x$interim$benefit
  ben$decision <- NA

  ben$decision[which(x$simresults$sel.dose == dose_1)] <- paste("Select", dose_1)
  ben$decision[which(x$simresults$sel.dose == dose_2)] <- paste("Select", dose_2)
  ben$decision[which(x$simresults$fut.trig)] <- "Futility"

  x_label <- paste("Estimated benefit", dose_1)
  y_label <- paste("Estimated benefit", dose_2)

  # Plot the bar chart with percentage labels
  ggplot(ben, aes(x = ben[,1], y = ben[,2], color = decision)) + geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    xlab(x_label) +
    ylab(y_label) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title = NULL))
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
