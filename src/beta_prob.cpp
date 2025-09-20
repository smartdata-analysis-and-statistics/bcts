#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double beta_prob_gt(double a, double b, double M, int n_draws) {
  NumericVector samples = Rcpp::rbeta(n_draws, a, b);
  return mean(samples > M);
}

//' @title Compute p-values for a t-distribution with Fixed Degrees of Freedom
//'
//' @description Simulates a single-arm binomial trial with a conjugate Beta prior,
//' and computes the Bayesian power - the proportion of trials where the posterior
//' probability that the response rate exceeds a threshold \code{M} is greater than or equal
//' to a specified decision threshold.
//'
//' @param B Integer. Number of trial simulations.
//' @param p_t Numeric in \[0, 1\]. True response probability for the treatment arm.
//' @param n_t Integer. Sample size of the treatment arm.
//' @param M Numeric in \[0, 1\]. Decision threshold on the response rate, e.g., \code{M = 0.6}.
//' @param threshold Numeric in \[0, 1\]. Posterior probability cutoff for declaring success,
//' e.g., \code{0.95}.
//' @param prior Character string. Either \code{"flat"} for a Beta(1,1) prior or
//' \code{"beta"} to specify a custom prior using \code{a_base} and \code{b_base}.
//' @param a_base Numeric. Alpha parameter for the Beta prior (only used if \code{prior = "beta"}).
//' @param b_base Numeric. Beta parameter for the Beta prior (only used if \code{prior = "beta"}).
//' @param n_draws Integer. Number of posterior draws per trial.
//' @param show_progress Logical. If \code{TRUE}, prints a simple progress bar to console.
//'
//' @return A list with: estimate (power), mc_se, successes, B.
//'
//' @examples
//' singlearm_beta_power(
//'   B = 1000, p_t = 0.75, n_t = 35, M = 0.60,
//'   threshold = 0.95, prior = "flat", n_draws = 2000
//' )
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
List singlearm_beta_power(int B,
                          double p_t,
                          int n_t,
                          double M,
                          double threshold,
                          std::string prior = "flat",
                          double a_base = 1,
                          double b_base = 1,
                          int n_draws = 2000,
                          bool show_progress = true) {

  int successes = 0;

  for (int b = 0; b < B; ++b) {
    // Simulate y_t ~ Binom(n_t, p_t)
    int y_t = R::rbinom(n_t, p_t);

    // Posterior Beta parameters
    double a_t, b_t;
    if (prior == "flat") {
      a_t = 1 + y_t;
      b_t = 1 + (n_t - y_t);
    } else if (prior == "beta") {
      a_t = a_base + y_t;
      b_t = b_base + (n_t - y_t);
    } else {
      stop("Invalid prior specification: must be 'flat' or 'beta'");
    }

    // Posterior samples for theta_t
    NumericVector draws = Rcpp::rbeta(n_draws, a_t, b_t);
    double prob = mean(draws > M);

    if (prob >= threshold) {
      successes++;
    }

    if (show_progress && (b + 1) % (B / 10) == 0) {
      Rcpp::Rcout << "."; // simple dot progress
    }
  }

  double est_power = static_cast<double>(successes) / B;
  double mc_se = std::sqrt(est_power * (1 - est_power) / B);

  return List::create(
    Named("estimate") = est_power,
    Named("mc_se") = mc_se,
    Named("B") = B,
    Named("successes") = successes
  );
}

//' @title Estimate Type-I Error for Single-Arm Trial
//'
//' @description Simulates a single-arm binomial trial under the null hypothesis,
//' and computes the type-I error: the proportion of simulations where posterior
//' Pr(θ > M) ≥ threshold even though θ = p_null.
//'
//' @param B Integer. Number of simulations.
//' @param n_t Integer. Sample size of the treatment arm.
//' @param M Numeric. Decision threshold for θ (on probability scale, e.g., 0.6).
//' @param threshold Posterior probability threshold γ (e.g., 0.95).
//' @param prior "flat" or "beta".
//' @param a_base Alpha parameter for Beta prior (if prior = "beta").
//' @param b_base Beta parameter for Beta prior (if prior = "beta").
//' @param n_draws Number of posterior draws per trial.
//' @param show_progress Logical. Show progress in console?
//'
//' @return A list with \code{estimate} (type-I error), \code{mc_se}, \code{B}, and \code{rejections}.
//'
//' @examples
//' singlearm_beta_type1(
//'   B = 1000, n_t = 35, M = 0.6,
//'   threshold = 0.95, prior = "flat", n_draws = 2000
//' )
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
List singlearm_beta_type1(int B,
                          int n_t,
                          double M,
                          double threshold,
                          std::string prior = "flat",
                          double a_base = 1,
                          double b_base = 1,
                          int n_draws = 2000,
                          bool show_progress = true) {

  int rejections = 0;

  for (int b = 0; b < B; ++b) {
    // Simulate response under null hypothesis: p_null = M
    int y_t = R::rbinom(n_t, M);

    double a_t, b_t;
    if (prior == "flat") {
      a_t = 1 + y_t;
      b_t = 1 + (n_t - y_t);
    } else if (prior == "beta") {
      a_t = a_base + y_t;
      b_t = b_base + (n_t - y_t);
    } else {
      stop("Invalid prior specification: must be 'flat' or 'beta'");
    }

    NumericVector draws = Rcpp::rbeta(n_draws, a_t, b_t);
    double prob = mean(draws > M);

    if (prob >= threshold) {
      rejections++;
    }

    if (show_progress && (b + 1) % (B / 10) == 0) {
      Rcpp::Rcout << ".";
    }
  }

  double type1 = static_cast<double>(rejections) / B;
  double mc_se = std::sqrt(type1 * (1 - type1) / B);

  return List::create(
    Named("estimate") = type1,
    Named("mc_se") = mc_se,
    Named("B") = B,
    Named("rejections") = rejections
  );
}
