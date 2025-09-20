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
```

### Launch the Shiny App

After installation, launch the interactive Shiny app by running:

```r
bcts::run_bcts_app()
```

This will open a browser window where you can:

* Simulate randomized and single-arm Bayesian trials,
* Calibrate posterior decision thresholds,
* Estimate Type-I error and power under different priors,
* Explore visualizations of posterior distributions and calibration traces.
 
The app supports trial design using **Bayesian Beta–Binomial models** and is intended to assist both applied users and methodologists.



