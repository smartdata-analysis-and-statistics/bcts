// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;

//' Estimate Power for Bayesian RCT Using Beta-Binomial Conjugate Model (C++)
//'
//' Performs Monte Carlo simulation to estimate the statistical power of a
//' Bayesian randomized controlled trial (RCT) using conjugate Beta-Binomial
//' models for the control and treatment arms. The function supports flat priors
//' and power priors on the control arm, and allows for fully customizable
//' baseline Beta prior parameters for both arms.
//'
//' A trial is declared successful if the posterior probability
//' \eqn{Pr(\theta_t - \theta_c > M \mid \text{data}) \ge \gamma}.
//'
//' @param B Integer. Number of simulated trials.
//' @param p_c,p_t Numeric in \[0, 1\]. True response probabilities in the control and treatment arms.
//' @param n_c,n_t Integers. Sample sizes in the control and treatment arms.
//' @param M Numeric. Margin on the risk-difference scale:
//'   - Negative for non-inferiority,
//'   - Zero for equivalence,
//'   - Positive for superiority.
//' @param threshold Numeric in (0, 1). Posterior probability threshold \eqn{\gamma}.
//' @param prior Character. Type of prior to use:
//'   - `"flat"`: independent Beta priors for both arms,
//'   - `"power"`: power prior on the control arm, flat prior on treatment.
//' @param prior_args List of prior hyperparameters. The following elements are supported:
//'   - `a0`: Discount factor for historical control data (only used if `prior = "power"`).
//'   - `y_0`, `n_0`: Number of responses and total patients in the historical control data.
//'   - `a_base_c`, `b_base_c`: Shape parameters of the Beta prior for the control arm.
//'   - `a_base_t`, `b_base_t`: Shape parameters of the Beta prior for the treatment arm.
//'     - Defaults for all `a_base_*` and `b_base_*` values are 1 (i.e., flat prior).
//' @param n_draws Integer. Number of posterior draws per trial for estimating the probability.
//' @param show_progress Logical. Show progress bar in the console.
//'
//' @return A logical vector of length `B`, indicating for each trial whether the decision
//' criterion was met (i.e., trial declared successful).
//'
//' @details This function uses Rcpp and vectorized binomial simulation to increase speed.
//' Posterior samples are drawn from Beta distributions parameterized using either flat priors
//' or power priors (for control) combined with observed trial data.
//'
//' @examples
//' prior_args <- list(
//'   a0 = 0.5, y_0 = 20, n_0 = 30,
//'   a_base_c = 1, b_base_c = 1,
//'   a_base_t = 2, b_base_t = 2
//' )
//' decisions <- rct_power_beta_binom_cpp_vec(
//'   B = 1000, p_c = 0.8, p_t = 0.8,
//'   n_c = 25, n_t = 25,
//'   M = -0.1, threshold = 0.9,
//'   prior = "power", prior_args = prior_args,
//'   n_draws = 2000, show_progress = FALSE
//' )
//' mean(decisions)  # Estimated power
//'
//' @export
// [[Rcpp::export]]
LogicalVector rct_power_beta_binom_cpp_vec(int B,
                                           double p_c, double p_t,
                                           int n_c, int n_t,
                                           double M,
                                           double threshold,
                                           std::string prior,
                                           List prior_args,
                                           int n_draws,
                                           bool show_progress = true) {
  LogicalVector decision(B);
  Progress prog(B, show_progress);

  IntegerVector y_c(B), y_t(B);
  for (int i = 0; i < B; ++i) {
    y_c[i] = R::rbinom(n_c, p_c);
    y_t[i] = R::rbinom(n_t, p_t);
  }

  // Prior hyperparameters
  double a0        = prior_args.containsElementNamed("a0")      ? as<double>(prior_args["a0"]) : 1.0;
  int y0           = prior_args.containsElementNamed("y_0")     ? as<int>(prior_args["y_0"])    : 0;
  int n0           = prior_args.containsElementNamed("n_0")     ? as<int>(prior_args["n_0"])    : 0;

  // Control arm baseline
  double a_base_c  = prior_args.containsElementNamed("a_base_c") ? as<double>(prior_args["a_base_c"]) : 1.0;
  double b_base_c  = prior_args.containsElementNamed("b_base_c") ? as<double>(prior_args["b_base_c"]) : 1.0;

  // Treatment arm baseline
  double a_base_t  = prior_args.containsElementNamed("a_base_t") ? as<double>(prior_args["a_base_t"]) : 1.0;
  double b_base_t  = prior_args.containsElementNamed("b_base_t") ? as<double>(prior_args["b_base_t"]) : 1.0;

  for (int b = 0; b < B; ++b) {
    if (Progress::check_abort()) break;

    double a_c, b_c, a_t, b_t;

    if (prior == "flat") {
      // Use arm-specific base values
      a_c = a_base_c + y_c[b];
      b_c = b_base_c + n_c - y_c[b];
      a_t = a_base_t + y_t[b];
      b_t = b_base_t + n_t - y_t[b];
    } else if (prior == "power") {
      // Control arm uses power prior
      a_c = a_base_c + a0 * y0 + y_c[b];
      b_c = b_base_c + a0 * (n0 - y0) + (n_c - y_c[b]);
      // Treatment arm uses flat prior with adjustable base
      a_t = a_base_t + y_t[b];
      b_t = b_base_t + n_t - y_t[b];
    } else {
      stop("Unknown prior: must be 'flat' or 'power'");
    }

    NumericVector theta_c = Rcpp::rbeta(n_draws, a_c, b_c);
    NumericVector theta_t = Rcpp::rbeta(n_draws, a_t, b_t);

    int count = sum((theta_t - theta_c) > M);
    decision[b] = (static_cast<double>(count) / n_draws) >= threshold;

    prog.increment();
  }

  return decision;
}
