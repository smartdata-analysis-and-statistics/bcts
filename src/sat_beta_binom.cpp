#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double beta_prob_gt(double a, double b, double M, int n_draws) {
  NumericVector samples = Rcpp::rbeta(n_draws, a, b);
  return mean(samples > M);
}


//' @title Analytic Power Calculation for Single-Arm Beta-Binomial Design
//'
//' @description Computes the exact Bayesian power for a single-arm binomial trial with a
//' conjugate Beta prior. Power is defined as the probability that the posterior probability
//' that the true response rate exceeds a threshold \code{M} is greater than or equal to
//' a prespecified cutoff \code{threshold}, under a fixed true response rate \code{p_t}.
//'
//' This function avoids simulation and uses a deterministic summation over all possible
//' outcomes of the binomial distribution.
//'
//' @param p_t Numeric in \[0, 1\]. True response probability for the treatment arm.
//' @param n_t Integer. Sample size of the treatment arm.
//' @param M Numeric in \[0, 1\]. Threshold on the response rate for decision-making,
//' e.g., \code{M = 0.6}.
//' @param threshold Numeric in \[0, 1\]. Posterior probability cutoff for declaring success,
//' e.g., \code{0.95}.
//' @param prior Character string specifying the prior distribution.
//' Options are:
//' "flat" for a non-informative Beta(1,1) prior;
//' "jeffreys" for the Jeffreys prior (Beta(0.5, 0.5));
//' "beta" for a custom Beta(\code{a_base}, \code{b_base}) prior (user must provide \code{a_base} and \code{b_base}).
//' @param a_base Numeric. Alpha parameter for the Beta prior (used only if \code{prior = "beta"}).
//' @param b_base Numeric. Beta parameter for the Beta prior (used only if \code{prior = "beta"}).
//'
//' @return A named list with:
//' \describe{
//'   \item{\code{estimate}}{The exact Bayesian power (a number between 0 and 1).}
//'   \item{\code{mc_se}}{\code{NA_real_}, returned for compatibility (Monte Carlo SE not applicable).}
//'   \item{\code{B}}{\code{NA_integer_}, returned for compatibility with simulation version.}
//'   \item{\code{successes}}{\code{NA_integer_}, returned for compatibility.}
//' }
//'
//' @examples
//' sat_betabinom_power_exact(
//'   p_t = 0.75, n_t = 35, M = 0.6,
//'   threshold = 0.95, prior = "flat"
//' )
//'
//' @seealso \code{\link{sat_betabinom_power}} for the simulation-based version.
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
List sat_betabinom_power_exact(double p_t,
                                int n_t,
                                double M,
                                double threshold,
                                std::string prior = "flat",
                                double a_base = 1,
                                double b_base = 1) {

  // Set default prior if flat
  if (prior == "flat") {
    a_base = 1.0;
    b_base = 1.0;
  } else if (prior == "jeffreys") {
    a_base = 0.5;
    b_base = 0.5;
  } else if (prior != "beta") {
    stop("Invalid prior specification: must be 'flat', 'jeffreys', or 'beta'");
  }

  double power = 0.0;

  for (int y = 0; y <= n_t; ++y) {
    double a_post = a_base + y;
    double b_post = b_base + n_t - y;

    // Posterior probability Pr(theta > M)
    double post_prob = 1.0 - R::pbeta(M, a_post, b_post, /*lower_tail=*/1, /*log_p=*/0);

    if (post_prob >= threshold) {
      // Pr(y | true p_t)
      double binom_prob = R::dbinom(y, n_t, p_t, /*give_log=*/0);
      power += binom_prob;
    }
  }

  return List::create(
    Named("estimate") = power,
    Named("mc_se") = R_NaReal,  // Not applicable for analytic
    Named("B") = NA_INTEGER,
    Named("successes") = NA_INTEGER
  );
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
//' @param prior Character string specifying the prior distribution.
//' Options are:
//' "flat" for a non-informative Beta(1,1) prior;
//' "jeffreys" for the Jeffreys prior (Beta(0.5, 0.5));
//' "beta" for a custom Beta(\code{a_base}, \code{b_base}) prior (user must provide \code{a_base} and \code{b_base}).
//' @param a_base Numeric. Alpha parameter for the Beta prior (only used if \code{prior = "beta"}).
//' @param b_base Numeric. Beta parameter for the Beta prior (only used if \code{prior = "beta"}).
//' @param show_progress Logical. If \code{TRUE}, prints a simple progress bar to console.
//'
//' @return A list with: estimate (power), mc_se, successes, B.
//'
//' @examples
//' sat_betabinom_power(
//'   B = 1000, p_t = 0.75, n_t = 35, M = 0.60,
//'   threshold = 0.95, prior = "flat"
//' )
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
List sat_betabinom_power(int B,
                          double p_t,
                          int n_t,
                          double M,
                          double threshold,
                          std::string prior = "flat",
                          double a_base = 1,
                          double b_base = 1,
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
    } else if (prior == "jeffreys") {
      a_t = 0.5 + y_t;
      b_t = 0.5 + (n_t - y_t);
    } else if (prior == "beta") {
      a_t = a_base + y_t;
      b_t = b_base + (n_t - y_t);
    } else {
      stop("Invalid prior specification: must be 'flat' or 'beta'");
    }

    // Posterior samples for theta_t
    //NumericVector draws = Rcpp::rbeta(n_draws, a_t, b_t);
    //double prob = mean(draws > M);
    double prob = 1.0 - R::pbeta(M, a_t, b_t, /*lower_tail=*/1, /*log_p=*/0);

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
//' @param p_null True response probability under H₀. Used to simulate Binomial outcomes.
//' @param M Numeric. Decision threshold for θ (on probability scale, e.g., 0.6).
//' @param threshold Posterior probability threshold γ (e.g., 0.95).
//' @param a_base Alpha parameter for Beta prior.
//' @param b_base Beta parameter for Beta prior.
//' @param show_progress Logical. Show progress in console?
//'
//' @return A list with \code{estimate} (type-I error), \code{mc_se}, \code{B}, and \code{rejections}.
//'
//' @examples
//' sat_betabinom_type1(
//'   B = 1000, n_t = 35, M = 0.6,
//'   threshold = 0.95, prior = "flat"
//' )
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
List sat_betabinom_type1(int B,
                         int n_t,
                         double p_null,        // True prob under null (simulate y ~ Bin(n_t, p_null))
                         double M,             // Decision margin (used in Pr(θ > M))
                         double threshold,     // Decision threshold for Pr(θ > M)
                         double a_base = 1,
                         double b_base = 1,
                         bool show_progress = true) {

  int rejections = 0;

  for (int b = 0; b < B; ++b) {
    // Simulate data under H₀: θ = p_null
    int y = R::rbinom(n_t, p_null);

    // Compute posterior shape
    double a_post = a_base + y;
    double b_post = b_base + (n_t - y);

    // Posterior tail probability Pr(θ > M)
    double pr_theta_gt_M = 1.0 - R::pbeta(M, a_post, b_post, 1, 0);

    // Check if posterior tail exceeds threshold (γ)
    if (pr_theta_gt_M >= threshold) {
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


//' @title Exact Type-I Error for Single-Arm Trial (Beta-Binomial)
//'
//' @description Computes the exact Type-I error for a single-arm binomial trial
//' using a conjugate Beta prior. This is done by summing the probability of all
//' outcomes where the posterior probability that \eqn{\theta > M} exceeds a
//' decision threshold \code{threshold}, under the null hypothesis.
//'
//' By default, the null hypothesis is \eqn{\theta = M}, but a more conservative
//' frequentist setting can be evaluated by setting \code{p_null < M}.
//'
//' @param n_t Integer. Sample size of the treatment arm.
//' @param M Numeric in \[0, 1\]. Decision threshold on the response rate, e.g., \code{M = 0.6}.
//' @param threshold Numeric in \[0, 1\]. Posterior probability cutoff for declaring success,
//' e.g., \code{threshold = 0.95}.
//' @param a_base Numeric. Alpha parameter for the Beta prior.
//' @param b_base Numeric. Beta parameter for the Beta prior.
//' @param p_null Optional. True response probability under the null hypothesis (e.g., \code{p_null = 0.6}).
//' If not specified, defaults to \code{p_null = M} (boundary case).
//'
//' @return A named list with:
//' \describe{
//'   \item{\code{estimate}}{Exact Type-I error (a number between 0 and 1).}
//'   \item{\code{mc_se}}{\code{NA_real_}, included for compatibility.}
//'   \item{\code{B}}{\code{NA_real_}, included for compatibility.}
//'   \item{\code{rejections}}{\code{NA_integer_}, included for compatibility.}
//' }
//'
//' @examples
//' sat_betabinom_type1_exact(n_t = 40, M = 0.65, threshold = 0.9)
//'
//' @seealso \code{\link{sat_betabinom_type1}} for the simulation-based version.
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
List sat_betabinom_type1_exact(int n_t,
                                double M,
                                double threshold,
                                double a_base = 1,
                                double b_base = 1,
                                double p_null = -1.0  // default to M if not provided
) {

  // Default: p_null = M if not specified
  if (p_null < 0.0) {
    p_null = M;
  }

  double type1 = 0.0;

  for (int y = 0; y <= n_t; ++y) {
    // Posterior parameters
    double a_post = a_base + y;
    double b_post = b_base + n_t - y;

    // Posterior probability Pr(θ > M)
    double post_prob = 1.0 - R::pbeta(M, a_post, b_post, /*lower_tail=*/1, /*log_p=*/0);

    // If threshold crossed, add to type I error (weighted by binomial Pr(y | p_null))
    if (post_prob >= threshold) {
      double binom_prob = R::dbinom(y, n_t, p_null, /*log=*/0);
      type1 += binom_prob;
    }
  }

  return List::create(
    Named("estimate") = type1,
    Named("mc_se") = R_NaReal,
    Named("B") = R_NaReal,
    Named("rejections") = NA_INTEGER
  );
}
