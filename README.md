<!-- badges: start -->

[![R-CMD-check](https://github.com/smartdata-analysis-and-statistics/bcts/actions/workflows/r.yml/badge.svg)](https://github.com/smartdata-analysis-and-statistics/bcts/actions/workflows/r.yml)
<!-- badges: end -->

## Overview

The **bcts** package provides functions to simulate, calibrate, and evaluate
**Bayesian non-inferiority (NI)** and **superiority trials** with binary outcomes.
It implements conjugate **Beta–Binomial models**, with support for:

- **Flat priors** (no external evidence),
- **Power priors** incorporating external control data.

The package is designed for both **methodological research** and **practical trial planning**,
offering tools to:

- Simulate trial operating characteristics (Type-I error, power, posterior probabilities),
- Calibrate posterior thresholds (γ) to target frequentist error rates (α),
- Evaluate borrowing strength from external data,
- Visualize calibration traces, prior/posterior distributions, and effective sample sizes.

A companion **Shiny app** is available to explore designs interactively.

---

## Installation

### Development version (from GitHub)

```r
# install.packages("remotes")
remotes::install_github("smartdata-analysis-and-statistics/bcts")
