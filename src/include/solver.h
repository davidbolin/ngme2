#ifndef __Solver__Solver__
#define __Solver__Solver__

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

#include "MatrixAlgebra.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cholmod.h>
#include <iostream>
#include <stdexcept>
#include <string.h>

#include <Eigen/CholmodSupport>
#ifdef USEMKL
#include <Eigen/PardisoSupport>
#endif

class sparse_llt_solver {
private:
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int>> R_eigen;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double, 0, int>> R_ldlt;
  Eigen::CholmodDecomposition<Eigen::SparseMatrix<double, 0, int>>
      R_cholmod_ldlt;
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double, 0, int>> R_supernodal;
#ifdef USEMKL
  Eigen::PardisoLLT<Eigen::SparseMatrix<double, 0, int>> R_pardiso;
  Eigen::PardisoLDLT<Eigen::SparseMatrix<double, 0, int>> R_pardiso_ldlt;
#endif
  int solver_type{0};
  int n{0}; // dimension of the factorized system (rows of Q)
  int N_iter{10};
  bool isSymmetric{true};
  Eigen::SparseMatrix<double, 0, int> Qi;
  Eigen::MatrixXd U, QU;
  bool Qi_computed{false}, QU_computed{false};
  // For non-symmetric mode we keep the last K to build normal equations and to
  // apply K^T on RHS when required
  Eigen::SparseMatrix<double, 0, int> K_last;

public:
  sparse_llt_solver() = default;
  sparse_llt_solver(int stype, int nin, int Ntrace, bool symmetric)
      : solver_type(stype), n(nin), N_iter(Ntrace), isSymmetric(symmetric) {
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  // Backward-compatible init plus symmetric toggle
  inline void init(int nin, int Ntrace, int /*max_iter*/, double /*tol*/,
                   int stype) {
    init(nin, Ntrace, /*symmetric*/ true, stype);
  }
  inline void init(int nin, int Ntrace, bool symmetric, int stype) {
    n = nin;
    N_iter = Ntrace;
    solver_type = stype;
    isSymmetric = symmetric;
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
    Qi_computed = QU_computed = false;
  }
  void analyze(const Eigen::SparseMatrix<double, 0, int> &M) {
    switch (solver_type) {
    case 0:
      if (isSymmetric) {
        R_eigen.analyzePattern(M);
      } else {
        R_eigen.analyzePattern(M.transpose() * M);
      }
      break;
    case 1:
      if (isSymmetric) {
        R_ldlt.analyzePattern(M);
      } else {
        R_ldlt.analyzePattern(M.transpose() * M);
      }
      break;
    case 2:
      if (isSymmetric) {
        R_supernodal.analyzePattern(M);
      } else {
        R_supernodal.analyzePattern(M.transpose() * M);
      }
      break;
    case 3:
      if (isSymmetric) {
        R_cholmod_ldlt.analyzePattern(M);
      } else {
        R_cholmod_ldlt.analyzePattern(M.transpose() * M);
      }
      break;
#ifdef USEMKL
    case 6:
      if (isSymmetric) {
        R_pardiso.analyzePattern(M);
      } else {
        R_pardiso.analyzePattern(M.transpose() * M);
      }
      break;
    case 7:
      if (isSymmetric) {
        R_pardiso_ldlt.analyzePattern(M);
      } else {
        R_pardiso_ldlt.analyzePattern(M.transpose() * M);
      }
      break;
#endif
    default:
      throw std::runtime_error("Pardiso solver not available (recompile with "
                               "USEMKL) or invalid solver_type");
    }
    QU_computed = false;
    n = isSymmetric ? M.rows() : M.cols();
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  void compute(const Eigen::SparseMatrix<double, 0, int> &M) {
    switch (solver_type) {
    case 0:
      if (isSymmetric) {
        R_eigen.factorize(M);
      } else {
        K_last = M;
        R_eigen.factorize(M.transpose() * M);
      }
      break;
    case 1:
      if (isSymmetric) {
        R_ldlt.factorize(M);
      } else {
        K_last = M;
        R_ldlt.factorize(M.transpose() * M);
      }
      break;
    case 2:
      if (isSymmetric) {
        R_supernodal.factorize(M);
      } else {
        K_last = M;
        R_supernodal.factorize(M.transpose() * M);
      }
      break;
    case 3:
      if (isSymmetric) {
        R_cholmod_ldlt.setMode(Eigen::CholmodLDLt);
        R_cholmod_ldlt.factorize(M);
      } else {
        K_last = M;
        R_cholmod_ldlt.setMode(Eigen::CholmodLDLt);
        R_cholmod_ldlt.factorize(M.transpose() * M);
      }
      break;
#ifdef USEMKL
    case 6:
      if (isSymmetric) {
        R_pardiso.factorize(M);
      } else {
        K_last = M;
        R_pardiso.factorize(M.transpose() * M);
      }
      break;
    case 7:
      if (isSymmetric) {
        R_pardiso_ldlt.factorize(M);
      } else {
        K_last = M;
        R_pardiso_ldlt.factorize(M.transpose() * M);
      }
      break;
#endif
    default:
      throw std::runtime_error("Pardiso solver not available (recompile with "
                               "USEMKL) or invalid solver_type");
    }
    Qi_computed = false;
    QU_computed = false;
    n = isSymmetric ? M.rows() : M.cols();
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  inline Eigen::ComputationInfo factorization_info() const {
    switch (solver_type) {
    case 0:
      return R_eigen.info();
    case 1:
      return R_ldlt.info();
    case 2:
      return R_supernodal.info();
    case 3:
      return R_cholmod_ldlt.info();
#ifdef USEMKL
    case 6:
      return R_pardiso.info();
    case 7:
      return R_pardiso_ldlt.info();
#endif
    default:
      return Eigen::InvalidInput;
    }
  }

  inline bool factorization_success() const {
    return factorization_info() == Eigen::Success;
  }

  // sample from N(Q^-1 mu, Q^-1), Q = G^T G + H^T H
  inline Eigen::VectorXd rMVN(SparseMatrix<double, 0, int> &G,
                              SparseMatrix<double, 0, int> &H,
                              Eigen::VectorXd &mu, Eigen::VectorXd &z1,
                              Eigen::VectorXd &z2) {
    Eigen::VectorXd x = G.transpose() * z1 + H.transpose() * z2 + mu;
    return solve(x);
  }

  inline Eigen::VectorXd solve(Eigen::VectorXd &v) {
    if (!isSymmetric && K_last.rows() > 0) {
      Eigen::VectorXd rhs = K_last.transpose() * v; // solve (K^T K) y = K^T v
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::VectorXd solve_raw(Eigen::VectorXd &rhs) {
    switch (solver_type) {
    case 0:
      return R_eigen.solve(rhs);
    case 1:
      return R_ldlt.solve(rhs);
    case 2:
      return R_supernodal.solve(rhs);
    case 3:
      return R_cholmod_ldlt.solve(rhs);
#ifdef USEMKL
    case 6:
      return R_pardiso.solve(rhs);
    case 7:
      return R_pardiso_ldlt.solve(rhs);
#endif
    default:
      throw std::runtime_error("Pardiso solver not available (recompile with "
                               "USEMKL) or invalid solver_type");
    }
  }

  inline Eigen::MatrixXd solve(Eigen::MatrixXd &v) {
    if (!isSymmetric && K_last.rows() > 0) {
      Eigen::MatrixXd rhs = K_last.transpose() * v;
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::MatrixXd solve_raw(Eigen::MatrixXd &rhs) {
    switch (solver_type) {
    case 0:
      return R_eigen.solve(rhs);
    case 1:
      return R_ldlt.solve(rhs);
    case 2:
      return R_supernodal.solve(rhs);
    case 3:
      return R_cholmod_ldlt.solve(rhs);
#ifdef USEMKL
    case 6:
      return R_pardiso.solve(rhs);
    case 7:
      return R_pardiso_ldlt.solve(rhs);
#endif
    default:
      throw std::runtime_error("Pardiso solver not available (recompile with "
                               "USEMKL) or invalid solver_type");
    }
  }

  inline Eigen::SparseMatrix<double, 0, int>
  solve(const Eigen::SparseMatrix<double, 0, int> &v) {
    if (!isSymmetric && K_last.rows() > 0) {
      Eigen::SparseMatrix<double, 0, int> rhs = K_last.transpose() * v;
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::SparseMatrix<double, 0, int>
  solve_raw(const Eigen::SparseMatrix<double, 0, int> &rhs) {
    switch (solver_type) {
    case 0:
      return R_eigen.solve(rhs);
    case 1:
      return R_ldlt.solve(rhs);
    case 2:
      return R_supernodal.solve(rhs);
    case 3:
      return R_cholmod_ldlt.solve(rhs);
#ifdef USEMKL
    case 6:
      return R_pardiso.solve(rhs);
    case 7:
      return R_pardiso_ldlt.solve(rhs);
#endif
    default:
      throw std::runtime_error("Pardiso solver not available (recompile with "
                               "USEMKL) or invalid solver_type");
    }
  }

  // Hutchinson estimator for tr(B Q^{-1} A Q^{-1}).
  // For symmetric K (Q=K), this equals tr(K^{-1} A K^{-1} B) by cyclicity.
  // For non-symmetric K (Q=K^T K), we internally use A_eff = K^T A so that
  //   tr(B Q^{-1} A_eff Q^{-1}) = tr((K^T K)^{-1} K^T A (K^T K)^{-1} B)
  double trace2(const Eigen::SparseMatrix<double, 0, int> &A,
                const Eigen::SparseMatrix<double, 0, int> &B,
                unsigned int seed = 0);
  double trace(const Eigen::SparseMatrix<double, 0, int> &,
               unsigned int seed = 0);
  double logdet() {
    switch (solver_type) {
    case 0:
      return log(R_eigen.determinant());
    case 1:
      return R_ldlt.vectorD().array().log().sum();
    case 2:
      return R_supernodal.logDeterminant();
    case 3:
      return R_cholmod_ldlt.logDeterminant();
#ifdef USEMKL
    case 6:
      throw std::runtime_error("Pardiso solver logdet not implemented");
    case 7:
      throw std::runtime_error("Pardiso LDLT solver logdet not implemented");
#endif
    default:
      throw std::runtime_error("Pardiso solver not available (recompile with "
                               "USEMKL) or invalid solver_type");
    }
  }
};

#endif
