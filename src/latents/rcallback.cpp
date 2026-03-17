#include "../operator.h"

// Var1: VAR(1) bivariate operator with Cayley reparameterization.
//
// Parameters: theta_K = (p1, p2, p3, p4), all unconstrained.
//
// Cayley transform:
//   J = [[0, p1], [-p1, 0]]           (skew-symmetric)
//   L = [[p2, 0], [p3, p4]]           (lower-triangular)
//   R = L L^T + eps * I               (positive definite)
//   S = J - R                         (stable: S + S^T = -2R < 0)
//   A = (I + S)(I - S)^{-1}          (spectral radius < 1 guaranteed)
//
// Precision matrix:
//   K = M0 + a11*M11 + a22*M22 + a12*M12 + a21*M21
//
// Numerical dK/dtheta_K is computed automatically by update_all().
// This implementation is fully thread-safe (no R callbacks).

using SM = SparseMatrix<double, 0, int>;

RCallback::RCallback(const Rcpp::List &operator_list)
    : Operator(operator_list),
      eps_cayley(1e-5)
{
  Rcpp::List mats = Rcpp::as<Rcpp::List>(operator_list["matrices"]);
  M0  = Rcpp::as<SM>(mats[0]);
  M11 = Rcpp::as<SM>(mats[1]);
  M22 = Rcpp::as<SM>(mats[2]);
  M12 = Rcpp::as<SM>(mats[3]);
  M21 = Rcpp::as<SM>(mats[4]);
}

void RCallback::build_KZ(const VectorXd &theta_K) {
  double p1 = theta_K(0), p2 = theta_K(1),
         p3 = theta_K(2), p4 = theta_K(3);

  // Skew-symmetric rotation part
  Eigen::Matrix2d J;
  J << 0.0, p1, -p1, 0.0;

  // Lower-triangular Cholesky factor
  Eigen::Matrix2d L;
  L << p2, 0.0, p3, p4;

  // Positive-definite contraction: R = L L^T + eps*I
  Eigen::Matrix2d R = L * L.transpose();
  R(0, 0) += eps_cayley;
  R(1, 1) += eps_cayley;

  // Stable matrix: all eigenvalues have Re < 0
  Eigen::Matrix2d I2 = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d S  = J - R;

  // Cayley transform: A = (I+S)(I-S)^{-1}
  Eigen::Matrix2d ImS = I2 - S;
  Eigen::Matrix2d IpS = I2 + S;
  Eigen::Matrix2d A   = ImS.inverse() * IpS;

  double a11 = A(0, 0), a12 = A(0, 1);
  double a21 = A(1, 0), a22 = A(1, 1);

  K = M0 + a11 * M11 + a22 * M22 + a12 * M12 + a21 * M21;
  // Z stays as identity (set by base constructor)
}
