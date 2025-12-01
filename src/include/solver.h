#ifndef __Solver__Solver__
#define __Solver__Solver__

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

#include <iostream>
#include <string.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MatrixAlgebra.h"

#include <cholmod.h>
#include <Eigen/CholmodSupport>
#ifdef __APPLE__
#include <Eigen/AccelerateSupport>
#endif
#ifdef USEMKL
#include <Eigen/PardisoSupport>
#endif

class sparse_llt_solver {
private:
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int> > R_eigen;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double, 0, int> > R_ldlt;
  Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double, 0, int> > R_simplicial;
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double, 0, int> > R_supernodal;
#ifdef __APPLE__
  Eigen::AccelerateLLT<Eigen::SparseMatrix<double, 0, int> > R_accelerate;
#endif
#ifdef USEMKL
  Eigen::PardisoLLT<Eigen::SparseMatrix<double, 0, int> > R_pardiso;
#endif
  int solver_type {0};
  int n {0}; // dimension of the factorized system (rows of Q)
  int N_iter {10};
  bool isSymmetric {true};
  Eigen::SparseMatrix<double, 0, int> Qi;
  Eigen::MatrixXd U, QU;
  bool Qi_computed {false}, QU_computed {false}; 
  // For non-symmetric mode we keep the last K to build normal equations and to apply K^T on RHS when required
  Eigen::SparseMatrix<double, 0, int> K_last;
public:
  sparse_llt_solver() = default;
  sparse_llt_solver(int stype, int nin, int Ntrace, bool symmetric)
    : solver_type(stype), n(nin), N_iter(Ntrace), isSymmetric(symmetric) { U.resize(n, N_iter); QU.resize(n, N_iter); }

  // Backward-compatible init plus symmetric toggle
  inline void init(int nin, int Ntrace, int /*max_iter*/, double /*tol*/, int stype) {
    init(nin, Ntrace, /*symmetric*/ true, stype);
  }
  inline void init(int nin, int Ntrace, bool symmetric, int stype) {
    n = nin; N_iter = Ntrace; solver_type = stype; isSymmetric = symmetric;
    U.resize(n, N_iter); QU.resize(n, N_iter);
    Qi_computed = QU_computed = false;
  }
  void analyze(const Eigen::SparseMatrix<double, 0, int> & M) {
    switch (solver_type) {
      case 0:
        if (isSymmetric) { R_eigen.analyzePattern(M); }
        else { R_eigen.analyzePattern(M.transpose() * M); }
        break;
      case 5:
        if (isSymmetric) { R_ldlt.analyzePattern(M); }
        else { R_ldlt.analyzePattern(M.transpose() * M); }
        break;
      case 1:
        if (isSymmetric) { R_simplicial.analyzePattern(M); }
        else { R_simplicial.analyzePattern(M.transpose() * M); }
        break;
      case 2:
        if (isSymmetric) { R_supernodal.analyzePattern(M); }
        else { R_supernodal.analyzePattern(M.transpose() * M); }
        break;
#ifdef __APPLE__
      case 3:
        if (isSymmetric) { R_accelerate.analyzePattern(M); }
        else { R_accelerate.analyzePattern(M.transpose() * M); }
        break;
#endif
#ifdef USEMKL
      case 4:
        if (isSymmetric) { R_pardiso.analyzePattern(M); }
        else { R_pardiso.analyzePattern(M.transpose() * M); }
        break;
#endif
      default:
        throw;
    }
    QU_computed = false;
    n = isSymmetric ? M.rows() : M.cols();
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  void compute(const Eigen::SparseMatrix<double, 0, int> & M) {
    switch (solver_type) {
      case 0:
        if (isSymmetric) { R_eigen.factorize(M); }
        else { K_last = M; R_eigen.factorize(M.transpose() * M); }
        break;
      case 5:
        if (isSymmetric) { R_ldlt.factorize(M); }
        else { K_last = M; R_ldlt.factorize(M.transpose() * M); }
        break;
      case 1:
        if (isSymmetric) { R_simplicial.factorize(M); }
        else { K_last = M; R_simplicial.factorize(M.transpose() * M); }
        break;
      case 2:
        if (isSymmetric) { R_supernodal.factorize(M); }
        else { K_last = M; R_supernodal.factorize(M.transpose() * M); }
        break;
#ifdef __APPLE__
      case 3:
        if (isSymmetric) { R_accelerate.factorize(M); }
        else { K_last = M; R_accelerate.factorize(M.transpose() * M); }
        break;
#endif
#ifdef USEMKL
      case 4:
        if (isSymmetric) { R_pardiso.factorize(M); }
        else { K_last = M; R_pardiso.factorize(M.transpose() * M); }
        break;
#endif
      default:
        throw;
    }
    Qi_computed = false;
    QU_computed = false;
    n = isSymmetric ? M.rows() : M.cols();
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  inline Eigen::ComputationInfo factorization_info() const {
    switch (solver_type) {
      case 0: return R_eigen.info();
      case 5: return R_ldlt.info();
      case 1: return R_simplicial.info();
      case 2: return R_supernodal.info();
#ifdef __APPLE__
      case 3: return R_accelerate.info();
#endif
#ifdef USEMKL
      case 4: return R_pardiso.info();
#endif
      default: return Eigen::InvalidInput;
    }
  }

  inline bool factorization_success() const {
    return factorization_info() == Eigen::Success;
  }

  // sample from N(Q^-1 mu, Q^-1), Q = G^T G + H^T H
  inline Eigen::VectorXd rMVN(
    SparseMatrix<double, 0, int> & G,
    SparseMatrix<double, 0, int> & H,
    Eigen::VectorXd & mu, 
    Eigen::VectorXd & z1, 
    Eigen::VectorXd & z2
  ) {
    Eigen::VectorXd x = G.transpose() * z1 + H.transpose() * z2 + mu;
    return solve(x);
  }


  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { 
    if (!isSymmetric && K_last.rows()>0) {
      Eigen::VectorXd rhs = K_last.transpose() * v; // solve (K^T K) y = K^T v
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::VectorXd solve_raw(Eigen::VectorXd &rhs) {
    switch (solver_type) {
      case 0:
        return R_eigen.solve(rhs);
      case 5:
        return R_ldlt.solve(rhs);
      case 1:
        return R_simplicial.solve(rhs);
      case 2:
        return R_supernodal.solve(rhs);
#ifdef __APPLE__
      case 3:
        return R_accelerate.solve(rhs);
#endif
#ifdef USEMKL
      case 4:
        return R_pardiso.solve(rhs);
#endif
      default:
        throw;
    }
  }

  inline Eigen::MatrixXd solve(Eigen::MatrixXd &v) { 
    if (!isSymmetric && K_last.rows()>0) {
      Eigen::MatrixXd rhs = K_last.transpose() * v;
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::MatrixXd solve_raw(Eigen::MatrixXd &rhs) { 
    switch (solver_type) {
      case 0:
        return R_eigen.solve(rhs);
      case 5:
        return R_ldlt.solve(rhs);
      case 1:
        return R_simplicial.solve(rhs);
      case 2:
        return R_supernodal.solve(rhs);
#ifdef __APPLE__
      case 3:
        return R_accelerate.solve(rhs);
#endif
#ifdef USEMKL
      case 4:
        return R_pardiso.solve(rhs);
#endif
      default:
        throw;
    }
  }

  inline Eigen::SparseMatrix<double, 0, int> solve(const Eigen::SparseMatrix<double, 0, int> &v) { 
    if (!isSymmetric && K_last.rows()>0) {
      Eigen::SparseMatrix<double,0,int> rhs = K_last.transpose() * v;
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::SparseMatrix<double, 0, int> solve_raw(const Eigen::SparseMatrix<double, 0, int> &rhs) { 
    switch (solver_type) {
      case 0:
        return R_eigen.solve(rhs);
      case 5:
        return R_ldlt.solve(rhs);
      case 1:
        return R_simplicial.solve(rhs);
      case 2:
        return R_supernodal.solve(rhs);
#ifdef __APPLE__
      case 3:
        return R_accelerate.solve(rhs);
#endif
#ifdef USEMKL
      case 4:
        return R_pardiso.solve(rhs);
#endif
      default:
        throw;
    }
  }
  
  // Hutchinson estimator for tr(B Q^{-1} A Q^{-1}).
  // For symmetric K (Q=K), this equals tr(K^{-1} A K^{-1} B) by cyclicity.
  // For non-symmetric K (Q=K^T K), we internally use A_eff = K^T A so that
  //   tr(B Q^{-1} A_eff Q^{-1}) = tr((K^T K)^{-1} K^T A (K^T K)^{-1} B)
  double trace2(const Eigen::SparseMatrix<double, 0, int> &A,
                const Eigen::SparseMatrix<double, 0, int> &B,
                unsigned int seed = 0);
  double trace(const Eigen::SparseMatrix<double, 0, int> &, unsigned int seed = 0);
  double logdet() {
    switch (solver_type) {
      case 0:
        return log(R_eigen.determinant());
      case 5:
        return R_ldlt.vectorD().array().log().sum();
      case 1:
        return R_simplicial.logDeterminant();
      case 2:
        return R_supernodal.logDeterminant();
#ifdef __APPLE__
      case 3:
        throw;
#endif
#ifdef USEMKL
      case 4:
        throw;
#endif
      default:
        throw;
    }
  }
};

#endif
