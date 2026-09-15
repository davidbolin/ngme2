#ifndef NGME_NIG_STD_H
#define NGME_NIG_STD_H

// Standardised coordinates for stationary NIG noise (Cabral, Bolin and Rue,
// 2023). In the native coordinates (theta_mu = mu, theta_sigma = log sigma,
// theta_nu = log nu) the variance sigma^2 + mu^2/nu is shared by all three
// parameters, so the likelihood has a long curved ridge along it; here that
// ridge is the first axis. Only the optimiser's coordinates change: the
// objective, the priors and the reported estimates stay native.
//
//   mode 0  native
//   mode 1  t = (log sigma_marg, zeta, log eta),  eta = 1/nu, zeta = mu/sigma,
//           sigma_marg = sqrt(sigma^2 + mu^2/nu)
//   mode 2  as 1, orthogonalised: zeta* = zeta sqrt(eta), eta* = eta / xi^2,
//           xi = 1 + zeta*^2 - |zeta*| sqrt(1 + zeta*^2), making the kurtosis
//           invariant to skewness
//   mode 3  as 2 with zeta* carried as asinh(zeta*)

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace nig_std {

// xi from Cabral sec 2.2. Because zeta* = zeta sqrt(eta) gives zeta^2 eta =
// zeta*^2, xi depends on zeta* alone, and lies in [1/2, 1].
inline double xi(double zstar) {
  const double z2 = zstar * zstar;
  return 1.0 + z2 - std::fabs(zstar) * std::sqrt(1.0 + z2);
}

// native (theta_mu, theta_sigma, theta_nu) -> t
inline Eigen::VectorXd from_native(int mode, const Eigen::VectorXd &native) {
  const double mu_v = native(0);
  const double sigma_v = std::exp(native(1));
  const double eta = std::exp(-native(2));
  const double zeta = mu_v / sigma_v;
  Eigen::VectorXd t(3);
  t(0) = std::log(std::sqrt(sigma_v * sigma_v + mu_v * mu_v * eta));
  if (mode == 1) {
    t(1) = zeta;
    t(2) = std::log(eta);
  } else {
    const double zstar = zeta * std::sqrt(eta);
    const double x = xi(zstar);
    t(1) = (mode == 3) ? std::asinh(zstar) : zstar;
    t(2) = std::log(eta / (x * x));
  }
  return t;
}

// t -> native (theta_mu, theta_sigma, theta_nu)
inline Eigen::VectorXd to_native(int mode, const Eigen::VectorXd &t) {
  const double sm = std::exp(t(0));
  double zeta, eta;
  if (mode == 1) {
    zeta = t(1);
    eta = std::exp(t(2));
  } else {
    // mode 3 puts zeta* on an asinh scale. xi -> 1/2 as |zeta*| grows, so the
    // excess kurtosis saturates at ~36 and the likelihood is flat out there;
    // asinh is linear near 0 and logarithmic in the tail, which matches that
    // sensitivity.
    const double zstar = (mode == 3) ? std::sinh(t(1)) : t(1);
    const double x = xi(zstar);
    eta = std::exp(t(2)) * x * x; // eta = eta* xi^2
    zeta = zstar / std::sqrt(eta);
  }
  const double D2 = 1.0 + zeta * zeta * eta;
  const double sigma_v = sm / std::sqrt(D2);
  Eigen::VectorXd out(3);
  out(0) = zeta * sigma_v;    // theta_mu = mu
  out(1) = std::log(sigma_v); // theta_sigma = log sigma
  out(2) = -std::log(eta);    // theta_nu = log nu = -log eta
  return out;
}

// d(native)/d(t) at the point with native coordinates `native`, by central
// differences on to_native. The map is a handful of flops, so this is far
// cheaper than the likelihood and avoids a hand-derivation for each mode.
inline Eigen::MatrixXd jacobian(int mode, const Eigen::VectorXd &native) {
  const Eigen::VectorXd t = from_native(mode, native);
  Eigen::MatrixXd J(3, 3);
  for (int j = 0; j < 3; ++j) {
    const double st = 1e-6 * std::max(1.0, std::fabs(t(j)));
    Eigen::VectorXd tp = t, tm = t;
    tp(j) += st;
    tm(j) -= st;
    J.col(j) = (to_native(mode, tp) - to_native(mode, tm)) / (2.0 * st);
  }
  return J;
}

} // namespace nig_std

#endif
