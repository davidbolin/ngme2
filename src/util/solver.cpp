#include "../include/solver.h"
#include <random>

using namespace Eigen;
double myround(double x) {
  if (x > 0) {
    return -1.0;
  } else {
    return 1.0;
  }
}

double sparse_llt_solver::trace(const SparseMatrix<double, 0, int> &M,
                                unsigned int seed) {
  if (QU_computed == 0) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < N_iter; ++j) {
        U(i, j) = dist(rng);
      }
    }
    U = U.unaryExpr(std::ref(myround));
    if (isSymmetric) {
      QU = solve(U);
    } else {
      // The LU path solves K x = u directly; the LLT/LDLT paths factorize
      // normal equations Q = K^T K and therefore need Q^{-1} u here.
      if (solver_type == 0)
        QU = R_eigen.solve(U);
      else if (solver_type == 1)
        QU = R_ldlt.solve(U);
      else if (solver_type == 2)
        QU = R_supernodal.solve(U);
      else if (solver_type == 3)
        QU = R_cholmod_ldlt.solve(U);
      else if (solver_type == 8)
        QU = R_lu.solve(U);
#ifdef __APPLE__
      else if (solver_type == 4)
        QU = R_accelerate.solve(U);
      else if (solver_type == 5)
        QU = R_accelerate_ldlt.solve(U);
#endif
#ifdef USEMKL
      else if (solver_type == 6)
        QU = R_pardiso.solve(U);
      else if (solver_type == 7)
        QU = R_pardiso_ldlt.solve(U);
#endif
      else
        throw std::runtime_error("Pardiso solver not available (recompile with "
                                 "USEMKL) or invalid solver_type");
    }
    QU_computed = 1;
  }

  Eigen::MatrixXd MQU;
  if (isSymmetric || K_last.rows() == 0 || solver_type == 8) {
    MQU = M * QU;
  } else {
    MQU = (K_last.transpose() * M) * QU;
  }
  double t = 0;
  for (int i = 0; i < N_iter; i++) {
    t += U.col(i).dot(MQU.col(i));
  }
  return t / N_iter;
}

// Hutchinson estimator for tr(K^{-1} A K^{-1} B).
// Symmetric case (Q=K):
//   Use QU = Q^{-1} U. Then S = Q^{-1}(A QU) via solve(A*QU).
//   Estimate tr = E[z^T B S z] = mean_i (u_i^T (B * S_i)).
// Non-symmetric case (Q=K^T K):
//   solve(X) internally applies K^T to the RHS before Q^{-1}, i.e.,
//   solve(X) = Q^{-1}(K^T X). Thus S = Q^{-1}(K^T A QU) which yields
//   tr(B Q^{-1} K^T A Q^{-1}) = tr((K^T K)^{-1} K^T A (K^T K)^{-1} B).
double sparse_llt_solver::trace2(const SparseMatrix<double, 0, int> &A,
                                 const SparseMatrix<double, 0, int> &B,
                                 unsigned int seed) {
  // Prepare Hutchinson vectors and their Q^{-1} images
  if (QU_computed == 0) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < N_iter; ++j) {
        U(i, j) = dist(rng);
      }
    }
    U = U.unaryExpr(std::ref(myround));
    if (isSymmetric) {
      QU = solve(U); // Q^{-1} U where Q=K
    } else {
      // The LU path solves K x = u directly; the LLT/LDLT paths factorize
      // normal equations Q = K^T K and therefore need Q^{-1} u here.
      if (solver_type == 0)
        QU = R_eigen.solve(U);
      else if (solver_type == 1)
        QU = R_ldlt.solve(U);
      else if (solver_type == 2)
        QU = R_supernodal.solve(U);
      else if (solver_type == 3)
        QU = R_cholmod_ldlt.solve(U);
      else if (solver_type == 8)
        QU = R_lu.solve(U);
#ifdef __APPLE__
      else if (solver_type == 4)
        QU = R_accelerate.solve(U);
      else if (solver_type == 5)
        QU = R_accelerate_ldlt.solve(U);
#endif
#ifdef USEMKL
      else if (solver_type == 6)
        QU = R_pardiso.solve(U);
      else if (solver_type == 7)
        QU = R_pardiso_ldlt.solve(U);
#endif
      else
        throw std::runtime_error("Pardiso solver not available (recompile with "
                                 "USEMKL) or invalid solver_type");
    }
    QU_computed = 1;
  }

  // First apply A to QU, then one more Q^{-1} using solve().
  // For LU on a non-symmetric K, solve() applies K^{-1} directly.
  // For the LLT/LDLT paths, solve() internally applies K^T to the RHS.
  Eigen::MatrixXd A_QU = A * QU;   // n x N_iter
  Eigen::MatrixXd S = solve(A_QU); // Q^{-1} A_eff Q^{-1} U
  Eigen::MatrixXd BS = B * S;      // n x N_iter

  double t = 0.0;
  for (int i = 0; i < N_iter; ++i) {
    t += U.col(i).dot(BS.col(i));
  }
  return t / static_cast<double>(N_iter);
}
