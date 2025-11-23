// fractional_operators.hpp
#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <string>

namespace rspde_cpp {

using SpMat = Eigen::SparseMatrix<double>;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

// Optional structure for rational-approximation roots
struct Roots {
  std::vector<double> rb;   // size m+1
  std::vector<double> rc;   // size m
  double factor;            // scalar scaling factor
};

// Provide your own implementation or link to a table-backed version.
Roots get_roots(int m, double beta, const std::string &type_interp = "linear");

// Core implementation that takes rb/rc/factor explicitly (for fractional beta)
std::pair<SpMat, SpMat> compute_fractional_operators_with_roots(
    const SpMat &C,
    const SpMat &Ci,
    const SpMat &G,
    double beta,
    int m,
    const VectorXd &tau,
    const VectorXd &theta_kappa,
    const MatrixXd &B_kappa,
    const std::vector<double> &rb,
    const std::vector<double> &rc,
    double roots_factor);

// Convenience wrapper that uses get_roots when beta is fractional.
std::pair<SpMat, SpMat> compute_fractional_operators(
    const SpMat &C,
    const SpMat &Ci,
    const SpMat &G,
    double beta,
    int m,
    const VectorXd &tau,
    const VectorXd &theta_kappa,
    const MatrixXd &B_kappa);

} // namespace rspde_cpp
