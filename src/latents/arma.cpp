#include "../operator.h"

// Build sparse shift matrix L (subdiagonal ones)
static Eigen::SparseMatrix<double,0,int> build_shift_L(int n) {
  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve(n-1);
  for (int i=1;i<n;i++) trip.emplace_back(i, i-1, 1.0);
  Eigen::SparseMatrix<double,0,int> L(n,n);
  L.setFromTriplets(trip.begin(), trip.end());
  return L;
}

// Build AR lag matrix C_j: -1 at (t, t-j) for t=j..n-1
static Eigen::SparseMatrix<double,0,int> build_Cj(int n, int j) {
  std::vector<Eigen::Triplet<double>> trip;
  for (int t=j; t<n; ++t) trip.emplace_back(t, t-j, -1.0);
  Eigen::SparseMatrix<double,0,int> C(n,n);
  C.setFromTriplets(trip.begin(), trip.end());
  return C;
}

Arma::Arma(const Rcpp::List& operator_list) :
  Operator(operator_list)
{
  // Dimension n from h
  n = Rcpp::as<Eigen::VectorXd>(operator_list["h"]).size();
  p = operator_list.containsElementNamed("p") ? Rcpp::as<int>(operator_list["p"]) : 1;
  q = operator_list.containsElementNamed("q") ? Rcpp::as<int>(operator_list["q"]) : 1;

  // Bases
  Cj.resize(std::max(0,p));
  for (int j=1; j<=p; ++j) Cj[j-1] = build_Cj(n, j);
  G.resize(n,n); G.setIdentity();
  L = build_shift_L(n);
  // Precompute powers of L up to max(p,q) for both AR(K) and MA(Z)
  int maxpq = std::max(p, q);
  Lpow.resize(std::max(0, maxpq));
  if (maxpq>0) {
    Lpow[0] = L;
    for (int k=2;k<=maxpq;k++) Lpow[k-1] = (Lpow[k-2] * L).pruned();
  }

  // Optional fixing masks
  // Only accept new names: fix_rho (AR) and fix_phi (MA)
  if (operator_list.containsElementNamed("fix_rho")) {
    Rcpp::LogicalVector fp = Rcpp::as<Rcpp::LogicalVector>(operator_list["fix_rho"]);
    fix_phi_mask.assign(fp.begin(), fp.end());
  }
  if (operator_list.containsElementNamed("fix_phi")) {
    Rcpp::LogicalVector ft = Rcpp::as<Rcpp::LogicalVector>(operator_list["fix_phi"]);
    fix_theta_mask.assign(ft.begin(), ft.end());
  }

  // Initialize K and Z from provided theta_K
  Eigen::VectorXd th = Rcpp::as<Eigen::VectorXd>(operator_list["theta_K"]);
  // If R passed actual coefficients, convert to raw (PACF space) here
  std::string param_space = operator_list.containsElementNamed("param_space") ?
      Rcpp::as<std::string>(operator_list["param_space"]) : std::string("raw");
  if (param_space == "coeff") {
    Eigen::VectorXd th_raw = th;
    auto clamp = [](double x){ return std::max(-0.999, std::min(0.999, x)); };
    auto atanh_safe = [&](double u){ u = clamp(u); return 0.5*std::log((1+u)/(1-u)); };
    // AR part (p<=2)
    if (p == 1 && th.size() >= 1) {
      double phi1 = clamp(th(0));
      th_raw(0) = 2.0 * atanh_safe(phi1); // pi1 = phi1
    } else if (p == 2 && th.size() >= 2) {
      double phi1 = clamp(th(0));
      double phi2 = clamp(th(1));
      double denom = std::max(1e-6, 1.0 - phi2);
      double pi2 = phi2;
      double pi1 = clamp(phi1 / denom);
      th_raw(0) = 2.0 * atanh_safe(pi1);
      th_raw(1) = 2.0 * atanh_safe(pi2);
    }
    // MA part (q<=2) lives after AR in theta_K
    if (q == 1 && th.size() >= p+1) {
      double theta1 = clamp(th(p+0));
      th_raw(p+0) = 2.0 * atanh_safe(theta1); // psi1 = theta1
    } else if (q == 2 && th.size() >= p+2) {
      double theta1 = clamp(th(p+0));
      double theta2 = clamp(th(p+1));
      double denom = std::max(1e-6, 1.0 - theta2);
      double psi2 = theta2;
      double psi1 = clamp(theta1 / denom);
      th_raw(p+0) = 2.0 * atanh_safe(psi1);
      th_raw(p+1) = 2.0 * atanh_safe(psi2);
    }
    th = th_raw;
  }
  // unified KZ build
  // AR part
  K.resize(n,n); K.setIdentity();
  if (p > 0) {
    double phi1 = 0.0, phi2 = 0.0;
    if (p == 1) {
      double pi1 = std::tanh(0.5 * th(0));
      phi1 = pi1;
    } else if (p == 2) {
      double pi1 = std::tanh(0.5 * th(0));
      double pi2 = std::tanh(0.5 * th(1));
      phi2 = pi2;
      phi1 = pi1 * (1.0 - pi2);
    }
    if (p >= 1) K += (-phi1) * Lpow[0];
    if (p >= 2) K += (-phi2) * Lpow[1];
  }
  // MA part
  Z.resize(n,n); Z.setIdentity();
  if (q > 0) {
    double th1 = 0.0, th2 = 0.0;
    if (q == 1) {
      th1 = std::tanh(0.5 * th(p + 0));
    } else if (q == 2) {
      double psi1 = std::tanh(0.5 * th(p + 0));
      double psi2 = std::tanh(0.5 * th(p + 1));
      th2 = psi2;
      th1 = psi1 * (1.0 - psi2);
    }
    if (q >= 1) Z += th1 * Lpow[0];
    if (q >= 2) Z += th2 * Lpow[1];
  }
}

// Implement base-required build_KZ using the same logic as in ctor
void Arma::build_KZ(const Eigen::VectorXd& theta_K) {
  // AR
  K.resize(n,n); K.setIdentity();
  if (p > 0) {
    double phi1 = 0.0, phi2 = 0.0;
    if (p == 1) {
      double pi1 = std::tanh(0.5 * theta_K(0));
      phi1 = pi1;
    } else if (p == 2) {
      double pi1 = std::tanh(0.5 * theta_K(0));
      double pi2 = std::tanh(0.5 * theta_K(1));
      phi2 = pi2;
      phi1 = pi1 * (1.0 - pi2);
    }
    if (p >= 1) K += (-phi1) * Lpow[0];
    if (p >= 2) K += (-phi2) * Lpow[1];
  }
  // Z
  Z.resize(n,n); Z.setIdentity();
  if (q > 0) {
    double th1 = 0.0, th2 = 0.0;
    if (q == 1) {
      th1 = std::tanh(0.5 * theta_K(p + 0));
    } else if (q == 2) {
      double psi1 = std::tanh(0.5 * theta_K(p + 0));
      double psi2 = std::tanh(0.5 * theta_K(p + 1));
      th2 = psi2;
      th1 = psi1 * (1.0 - psi2);
    }
    if (q >= 1) Z += th1 * Lpow[0];
    if (q >= 2) Z += th2 * Lpow[1];
  }
}
