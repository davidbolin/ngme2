#include "prior.h"

// log(pi(value; params))
double PriorUtil::log_dens(
  const string& type,
  const VectorXd& params,
  double value
) {
  if (type == "none") {
    return 0.0;
  } else if (type == "normal") {
    double mu = params(0);
    double prec = params(1);
    double sigma = 1.0 / std::sqrt(prec);
    return R::dnorm(value, mu, sigma, true);
  } else if (type == "inv.exponential") {
    double lambda = params(0);
    double lower = params.size() > 1 ? params(1) : 0.0;
    // value is theta on unconstrained scale, nu = lower + exp(theta)
    double exp_value = std::exp(value);
    double nu = lower + exp_value;
    return std::log(lambda) - lambda / nu - 2.0 * std::log(nu) + value;
  } else {
    throw std::invalid_argument("Unknown prior type");
  }
}

// d log(pi(value; params)) / d value
double PriorUtil::d_log_dens(
  const string& type,
  const VectorXd& params,
  double value
) {
  if (type == "none") {
    return 0;
  } else if (type == "normal") {
    double mu = params(0);
    double prec = params(1);
    double sigma = 1.0 / std::sqrt(prec);
    return (mu-value) / (sigma*sigma);
  } else if (type == "pc.sd") {
    // value = log(sd)
    // log(1/sd^2) = log(prec)
    // pi(log.prec) = pi(log.sd) * |1/2|
    double lambda = params(0);
    double tmp = exp(-value/2);
    return 2 * 0.25 * exp(-value - tmp*lambda) * lambda * (lambda - tmp);
  } else if (type == "half.cauchy") {
    // sigma ~ half-cauchy(0, gamma)
    // compute log_sigma
    double gamma = params(0);
    return (-exp(2*value) + gamma*gamma) / 
      (exp(2*value) + gamma*gamma);
  } else if (type == "inv.exponential") {
    double lambda = params(0);
    double lower = params.size() > 1 ? params(1) : 0.0;
    // value is theta, nu = lower + exp(theta)
    double exp_value = std::exp(value);
    double nu = lower + exp_value;
    return lambda * exp_value / (nu * nu) - 2.0 * exp_value / nu + 1.0;
  } else {
    throw std::invalid_argument("Unknown prior type");
  }
}
