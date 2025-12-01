/*
    Matern model (stationary or non-stationary kappa):
        kappa(s) = exp(B_theta_K * theta_K)
        alpha is the smoothness parameter
        K = kappa^2 * C + G for integer alpha (with optional fractional case)
*/

#include "../operator.h"
#include "fractional/fractional_operators.hpp"
#include <stdexcept>

Matern::Matern(const Rcpp::List &operator_list)
    : Operator(operator_list),
      G(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["G"])),
      C(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["C"])),
      Ci(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Ci"])),
      Bk_dense(Rcpp::as<MatrixXd>(operator_list["B_K"])),
      alpha(Rcpp::as<double>(operator_list["alpha"])), Cdiag(C.diagonal()) {
  // Stationarity flag and basis matrix for log kappa (dense)
  stationary = operator_list.containsElementNamed("stationary")
                   ? Rcpp::as<bool>(operator_list["stationary"])
                   : true;

  // Optional fractional controls
  fix_alpha = operator_list.containsElementNamed("fix_alpha")
                  ? Rcpp::as<bool>(operator_list["fix_alpha"])
                  : true;
  m = operator_list.containsElementNamed("rational_order")
          ? Rcpp::as<int>(operator_list["rational_order"])
          : 0;
  dim = operator_list.containsElementNamed("spatial_dim")
            ? Rcpp::as<int>(operator_list["spatial_dim"])
            : 2;

  const int n = G.rows();
  if (Bk_dense.rows() != n) {
    Rcpp::stop("B_K has %d rows but expected %d",
               static_cast<int>(Bk_dense.rows()), n);
  }
}

void Matern::build_KZ(const VectorXd &theta_K) {
  using namespace rspde_cpp;
  // theta_K layout:
  //   if (!fix_alpha): [eta_alpha, theta_K ...]
  //   if ( fix_alpha): [theta_K ...]
  int offset = 0;
  if (!fix_alpha) {
    double eta_alpha = theta_K(0);
    double L = 0.5 * static_cast<double>(dim);
    double sig = 1.0 / (1.0 + std::exp(-eta_alpha));
    alpha = L + (4.0 - L) * sig; // (L,4)
    offset = 1;
  }

  const int n_kappa = static_cast<int>(theta_K.size()) - offset;
  if (n_kappa <= 0) {
    throw std::invalid_argument("theta_K must contain kappa coefficients");
  }
  if (Bk_dense.cols() != n_kappa) {
    Rcpp::stop("Length of theta_K (%d) does not match B_K columns (%d)",
               n_kappa, static_cast<int>(Bk_dense.cols()));
  }
  VectorXd theta_kappa = theta_K.segment(offset, n_kappa);

  // kappa(s) = exp(B_theta_K * theta_kappa)
  VectorXd log_kappa = Bk_dense * theta_kappa;
  VectorXd kappa = log_kappa.array().exp();
  VectorXd kappa2 = kappa.array().square();

  if (std::abs(alpha - 2) < 1e-6) {
    // Integer case alpha = 2
    // K = G + C * diag(kappa^2)
    SparseMatrix<double> KCK = (C * kappa2.asDiagonal()).eval();
    K = (G + KCK);
    Z.setIdentity();
  } else if (std::abs(alpha - 4) < 1e-6) {
    // Integer case alpha = 4
    SparseMatrix<double> KCK = (C * kappa2.asDiagonal()).eval();
    K = (G + KCK) * Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
    Z.setIdentity();
  } else {
    // Fractional case
    double beta = alpha / 2.0;
    // tau defaults to 1 for now
    VectorXd tau(1);
    tau(0) = 1.0;
    auto pairKZ = compute_fractional_operators(C, Ci, G, beta, m, tau,
                                               theta_kappa, Bk_dense);
    K = pairKZ.first;  // Pl
    Z = pairKZ.second; // Pr
  }
}
