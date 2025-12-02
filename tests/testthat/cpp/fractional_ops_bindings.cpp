#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// We include the implementation directly to simplify testing setup.
// Using angle brackets to prevent sourceCpp from automatically compiling the
// file separately (which causes duplicate symbols), while still including the
// code.
#include <fractional_operators.cpp>

using namespace Rcpp;
using namespace rspde_cpp;

// [[Rcpp::export]]
List cpp_fractional_ops_bind(const Eigen::SparseMatrix<double> &C,
                             const Eigen::SparseMatrix<double> &Ci,
                             const Eigen::SparseMatrix<double> &G, double beta,
                             int m, const Eigen::VectorXd &tau,
                             const Eigen::VectorXd &theta_kappa,
                             const Eigen::SparseMatrix<double> &B_kappa) {

  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> res =
      compute_fractional_operators(C, Ci, G, beta, m, tau, theta_kappa,
                                   B_kappa);

  return List::create(Named("Pl") = res.first, Named("Pr") = res.second);
}

// [[Rcpp::export]]
List cpp_fractional_ops_with_roots_bind(
    const Eigen::SparseMatrix<double> &C, const Eigen::SparseMatrix<double> &Ci,
    const Eigen::SparseMatrix<double> &G, double beta, int m,
    const Eigen::VectorXd &tau, const Eigen::VectorXd &theta_kappa,
    const Eigen::SparseMatrix<double> &B_kappa, const std::vector<double> &rb,
    const std::vector<double> &rc, double roots_factor) {

  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> res =
      compute_fractional_operators_with_roots(
          C, Ci, G, beta, m, tau, theta_kappa, B_kappa, rb, rc, roots_factor);

  return List::create(Named("Pl") = res.first, Named("Pr") = res.second);
}
