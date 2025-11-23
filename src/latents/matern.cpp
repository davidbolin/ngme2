/*
    Matern model with stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = kappa
        K = kappa^2 * C + G
*/

#include "../operator.h"
#include "fractional/fractional_operators.hpp"

Matern::Matern(const Rcpp::List &operator_list)
    : Operator(operator_list),
      G(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["G"])),
      C(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["C"])),
      Ci(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Ci"])),
      alpha(Rcpp::as<double>(operator_list["alpha"])), Cdiag(C.diagonal()) {
  // Optional fractional controls
  fix_alpha = operator_list.containsElementNamed("fix_alpha")
                  ? Rcpp::as<bool>(operator_list["fix_alpha"])
                  : true;
  m = operator_list.containsElementNamed("m")
          ? Rcpp::as<int>(operator_list["m"])
          : 0;
  dim = operator_list.containsElementNamed("spatial_dim")
            ? Rcpp::as<int>(operator_list["spatial_dim"])
            : 2;
}

void Matern::build_KZ(const VectorXd &theta_K) {
  using namespace rspde_cpp;
  const int n = G.rows();
  // theta_K layout:
  //   if (!fix_alpha): [eta_alpha, log_kappa]
  //   if ( fix_alpha): [log_kappa]
  int offset = 0;
  if (!fix_alpha) {
    double eta_alpha = theta_K(0);
    double L = 0.5 * static_cast<double>(dim);
    double sig = 1.0 / (1.0 + std::exp(-eta_alpha));
    alpha = L + (4.0 - L) * sig; // (L,4)
    offset = 1;
  }
  Eigen::SparseMatrix<double> Bk(n, 1);
  Bk.reserve(n);
  for (int i = 0; i < n; ++i)
    Bk.insert(i, 0) = 1.0;
  Bk.makeCompressed();
  VectorXd theta(1);
  theta(0) = (theta_K.size() > offset) ? theta_K(offset) : 0.0; // log kappa
  // tau defaults to 1
  VectorXd tau(1);
  tau(0) = 1.0;

  bool is_integer = std::abs(alpha - std::round(alpha)) < 1e-6;
  if (is_integer) {
    double kappa = std::exp(theta(0));
    SparseMatrix<double> KCK = kappa * kappa * C;
    int ialpha = static_cast<int>(std::round(alpha));
    if (ialpha == 2) {
      K = (G + KCK);
    } else if (ialpha == 4) {
      K = (G + KCK) * Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
    } else {
      throw("alpha not equal to 2 or 4 is not implemented");
    }
    Z.setIdentity();
  } else {
    double beta = alpha / 2.0;
    std::pair<rspde_cpp::SpMat, rspde_cpp::SpMat> pairKZ;
    std::cout << " beta = " << beta << std::endl;
    std::cout << " theta = " << theta << std::endl;
    pairKZ = compute_fractional_operators(C, Ci, G, beta, 1, tau, theta, Bk);
    K = pairKZ.first;  // Pl
    Z = pairKZ.second; // Pr
  }
}
