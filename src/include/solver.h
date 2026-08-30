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
#include <Eigen/SparseLU>
#ifdef __APPLE__
#include <Eigen/AccelerateSupport>
#endif
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
#ifdef __APPLE__
  Eigen::AccelerateLLT<Eigen::SparseMatrix<double, 0, int>> R_accelerate;
  Eigen::AccelerateLDLT<Eigen::SparseMatrix<double, 0, int>> R_accelerate_ldlt;
#endif
#ifdef USEMKL
  Eigen::PardisoLLT<Eigen::SparseMatrix<double, 0, int>> R_pardiso;
  Eigen::PardisoLDLT<Eigen::SparseMatrix<double, 0, int>> R_pardiso_ldlt;
#endif
  // Direct LU of K, used instead of a Cholesky of K^T K when the operator is
  // not symmetric.
  Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int>>
      R_lu;
  bool use_lu{false};
  // Whether the last R_lu.factorize() succeeded. Eigen's SparseLU returns from
  // factorize() *before* it sets up m_Lstore when the matrix is singular, and
  // its solve path only guards that with an eigen_assert, which is compiled
  // out here (NDEBUG). Solving through such a decomposition dereferences
  // uninitialised pointers and crashes the R session, so every LU solve has to
  // check this flag first.
  bool lu_ok{false};

  void require_lu() const {
    if (!lu_ok)
      throw std::runtime_error(
          "LU factorization of the operator matrix K failed: K is singular "
          "or numerically rank-deficient");
  }

  int solver_type{0};
  int n{0}; // dimension of the factorized system (rows of Q)
  int N_iter{10};
  bool isSymmetric{true};
  Eigen::SparseMatrix<double, 0, int> Qi;
  Eigen::MatrixXd U, QU;
  bool Qi_computed{false}, QU_computed{false};
  // Hutchinson probe vectors.
  bool U_computed{false};
  unsigned int U_seed{0};
  // Set when U is the scaled identity rather than random probes, i.e. when the
  // probe budget is at least the dimension (n <= N_iter) and the trace is
  // therefore computed exactly.
  bool exact_trace{false};
  // Sample variance of the individual Hutchinson probe values from the most
  // recent trace()/trace2() call. The estimator averages N_iter probes, so the
  // variance it contributes to the gradient is last_probe_var_ / N_iter. This
  // is what lets the probe count be tuned against the Gibbs noise instead of
  // being fixed a priori. Zero when the trace was taken exactly.
  double last_probe_var_{0.0};
  // Selected (Takahashi) inverse of the factorized matrix, held on the
  // sparsity pattern of the Cholesky factor. Exact where it is defined, and
  // for a low-fill factor far cheaper than a Hutchinson estimate: it is built
  // once per factorization and then reused by every trace in that iteration.
  Eigen::SparseMatrix<double, 0, int> S_sel;
  bool S_sel_ready{false};
  // The default backends (CHOLMOD, Accelerate) do not expose their Cholesky
  // factor, so the selected inverse keeps its own simplicial factorization.
  // That is an extra factorization, but it is only ever taken when the factor
  // is low-fill, where it costs far less than the probe solves it replaces.
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int>> selinv_llt;
  Eigen::SparseMatrix<double, 0, int> M_sym_;
  bool M_sym_ready{false};
  bool selinv_llt_ready{false};
  int selinv_viable_{-1}; // -1 undecided, 0 fill too high, 1 usable
  void ensure_U(unsigned int seed);
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
  // nonsym_mode: 0 = LU of K (default), 1 = Cholesky of the normal equations
  // K^T K (the historical path). Chosen by control_opt(nonsym_solver = ).
  inline void init(int nin, int Ntrace, bool symmetric, int stype,
                   int nonsym_mode = 0) {
    use_lu = !symmetric && nonsym_mode == 0;
    lu_ok = false;
    n = nin;
    N_iter = Ntrace;
    solver_type = stype;
    isSymmetric = symmetric;
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
    Qi_computed = QU_computed = false;
    S_sel_ready = false;
    U_computed = false;
    exact_trace = false;
  }
  void analyze(const Eigen::SparseMatrix<double, 0, int> &M) {
    if (use_lu) {
      R_lu.analyzePattern(M);
      lu_ok = false;
      QU_computed = false;
    S_sel_ready = false;
      return;
    }
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
#ifdef __APPLE__
    case 4:
      if (isSymmetric) {
        R_accelerate.analyzePattern(M);
      } else {
        R_accelerate.analyzePattern(M.transpose() * M);
      }
      break;
    case 5:
      if (isSymmetric) {
        R_accelerate_ldlt.analyzePattern(M);
      } else {
        R_accelerate_ldlt.analyzePattern(M.transpose() * M);
      }
      break;
#endif
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
    S_sel_ready = false;
    n = isSymmetric ? M.rows() : M.cols();
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  void compute(const Eigen::SparseMatrix<double, 0, int> &M) {
    if (!use_lu && isSymmetric) {
      M_sym_ = M;
      M_sym_ready = true;
      selinv_llt_ready = false;
    }
    if (use_lu) {
      K_last = M;
      R_lu.factorize(M);
      lu_ok = (R_lu.info() == Eigen::Success);
      Qi_computed = false;
      QU_computed = false;
    S_sel_ready = false;
      const int new_n = M.cols();
      if (new_n != n)
        U_computed = false;
      n = new_n;
      U.resize(n, N_iter);
      QU.resize(n, N_iter);
      return;
    }
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
#ifdef __APPLE__
    case 4:
      if (isSymmetric) {
        R_accelerate.factorize(M);
      } else {
        K_last = M;
        R_accelerate.factorize(M.transpose() * M);
      }
      break;
    case 5:
      if (isSymmetric) {
        R_accelerate_ldlt.factorize(M);
      } else {
        K_last = M;
        R_accelerate_ldlt.factorize(M.transpose() * M);
      }
      break;
#endif
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
    S_sel_ready = false;
    const int new_n = isSymmetric ? M.rows() : M.cols();
    if (new_n != n)
      U_computed = false; // resize below invalidates the probe vectors
    n = new_n;
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
  }

  inline Eigen::ComputationInfo factorization_info() const {
    if (use_lu)
      return R_lu.info();
    switch (solver_type) {
    case 0:
      return R_eigen.info();
    case 1:
      return R_ldlt.info();
    case 2:
      return R_supernodal.info();
    case 3:
      return R_cholmod_ldlt.info();
#ifdef __APPLE__
    case 4:
      return R_accelerate.info();
    case 5:
      return R_accelerate_ldlt.info();
#endif
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

  // True when K is factorized by SparseLU rather than by a Cholesky. Only that
  // path is left in an unusable state by a failed factorization -- see lu_ok.
  inline bool uses_lu() const { return use_lu; }

  // sample from N(Q^-1 mu, Q^-1), Q = G^T G + H^T H
  inline Eigen::VectorXd rMVN(const SparseMatrix<double, 0, int> &G,
                              const SparseMatrix<double, 0, int> &H,
                              Eigen::VectorXd &mu, Eigen::VectorXd &z1,
                              Eigen::VectorXd &z2) {
    Eigen::VectorXd x = G.transpose() * z1 + H.transpose() * z2 + mu;
    return solve(x);
  }

  inline Eigen::VectorXd solve(Eigen::VectorXd &v) {
    if (use_lu) {
      require_lu();
      return R_lu.solve(v);
    }
    if (!isSymmetric && K_last.rows() > 0) {
      Eigen::VectorXd rhs = K_last.transpose() * v; // solve (K^T K) y = K^T v
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::VectorXd solve_raw(Eigen::VectorXd &rhs) {
    if (use_lu) {
      require_lu();
      return R_lu.solve(rhs);
    }
    switch (solver_type) {
    case 0:
      return R_eigen.solve(rhs);
    case 1:
      return R_ldlt.solve(rhs);
    case 2:
      return R_supernodal.solve(rhs);
    case 3:
      return R_cholmod_ldlt.solve(rhs);
#ifdef __APPLE__
    case 4:
      return R_accelerate.solve(rhs);
    case 5:
      return R_accelerate_ldlt.solve(rhs);
#endif
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
    if (use_lu) {
      require_lu();
      return R_lu.solve(v);
    }
    if (!isSymmetric && K_last.rows() > 0) {
      Eigen::MatrixXd rhs = K_last.transpose() * v;
      return solve_raw(rhs);
    }
    return solve_raw(v);
  }

  inline Eigen::MatrixXd solve_raw(Eigen::MatrixXd &rhs) {
    if (use_lu) {
      require_lu();
      return R_lu.solve(rhs);
    }
    switch (solver_type) {
    case 0:
      return R_eigen.solve(rhs);
    case 1:
      return R_ldlt.solve(rhs);
    case 2:
      return R_supernodal.solve(rhs);
    case 3:
      return R_cholmod_ldlt.solve(rhs);
#ifdef __APPLE__
    case 4:
      return R_accelerate.solve(rhs);
    case 5:
      return R_accelerate_ldlt.solve(rhs);
#endif
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
#ifdef __APPLE__
    case 4:
      return R_accelerate.solve(rhs);
    case 5:
      return R_accelerate_ldlt.solve(rhs);
#endif
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
  // Variance of the individual probe values in the last trace call (0 if exact).
  double last_probe_var() const { return last_probe_var_; }
  // Simplicial factorizations expose matrixL(); the CHOLMOD ones do not, and
  // the LU path is not a Cholesky at all.
  bool selinv_supported() const { return !use_lu && isSymmetric; }
  // nnz(L)/n. The cost of the selected inverse scales with the fill of the
  // factor, so this is what decides whether it beats probing.
  bool ensure_selinv_factor();
  double fill_ratio();
  bool build_selinv();
  // tr(Q^{-1} M) from the selected inverse. False when an entry M needs falls
  // outside the factor's pattern, leaving the caller to fall back to probing.
  bool selinv_trace(const Eigen::SparseMatrix<double, 0, int> &M, double &out);
  int get_N_iter() const { return N_iter; }
  bool is_exact_trace() const { return exact_trace; }
  // Change the probe budget at run time. Invalidates the cached probes so the
  // next trace call regenerates them (and switches to the exact path if the
  // budget now reaches the dimension).
  void set_N_iter(int Ntrace) {
    if (Ntrace < 1 || Ntrace == N_iter)
      return;
    N_iter = Ntrace;
    U.resize(n, N_iter);
    QU.resize(n, N_iter);
    U_computed = false;
    QU_computed = false;
    S_sel_ready = false;
    exact_trace = false;
  }

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
#ifdef __APPLE__
    case 4:
      throw std::runtime_error("Accelerate solver not available");
    case 5:
      throw std::runtime_error("Accelerate LDLT solver not available");
#endif
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
