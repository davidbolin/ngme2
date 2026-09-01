#ifndef __Solver__Solver__
#define __Solver__Solver__

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

#include "MatrixAlgebra.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cholmod.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string.h>

#include <Eigen/CholmodSupport>
#include <Eigen/SparseLU>
#ifdef __APPLE__
#include <Eigen/AccelerateSupport>
#include <unistd.h>
#endif
#ifdef USEMKL
#include <Eigen/PardisoSupport>
#endif

#ifdef __APPLE__
// ---------------------------------------------------------------------------
// Apple's Accelerate sparse solvers are NOT fork-safe.
//
// Eigen's AccelerateSupport wraps Apple's Sparse Solvers, which parallelise
// internally through libdispatch (GCD). GCD is documented as unusable in a
// process forked from one that has already initialised it: the child must
// exec() first. R's `parallel::mclapply()` forks and does not exec, so ANY use
// of an Accelerate factorization inside an mclapply worker aborts the child --
// and, because R's fork workers share the session's error handling, takes the
// whole R session down with it. Measured on this build: forked workers survive
// 3/3 with `solver_backend = "cholmod"` and 3/3 with `"eigen"`, and abort with
// `"accelerate"`, for Gaussian and non-Gaussian models alike.
//
// This is not something the caller can be expected to know, and the failure is
// a crash rather than an error, so the library declines to walk into it: a
// solver constructed in a forked child silently uses the equivalent CHOLMOD
// factorization instead. The two compute the same Cholesky; only the library
// differs, so results are unchanged up to floating-point associativity, and a
// forked worker is doing per-fold work where the solver's speed is not what
// dominates anyway.
//
// The child is detected by comparing the current pid against the one recorded
// when the shared library was loaded. `inline` gives one shared copy across
// translation units (C++17), initialised at load time.
inline const pid_t ngme_origin_pid = ::getpid();
inline bool ngme_in_forked_child() { return ::getpid() != ngme_origin_pid; }

// Map an Accelerate backend index to its CHOLMOD equivalent when -- and only
// when -- we are running in a forked child. 4 (Accelerate LLT) -> 2 (CHOLMOD
// supernodal LLT); 5 (Accelerate LDLT) -> 3 (CHOLMOD LDLT). Every other index
// is returned unchanged.
inline int ngme_fork_safe_stype(int stype) {
  if (!ngme_in_forked_child() || (stype != 4 && stype != 5))
    return stype;
  static thread_local bool warned = false;
  if (!warned) {
    warned = true;
    std::fprintf(stderr,
                 "ngme2: Apple Accelerate is not fork-safe; this forked worker "
                 "is using CHOLMOD instead. Use a PSOCK cluster "
                 "(parallel::makePSOCKcluster) to keep Accelerate.\n");
  }
  return (stype == 4) ? 2 : 3;
}
#else
inline int ngme_fork_safe_stype(int stype) { return stype; }
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
  // Held by pointer so disable_selinv() can hand the factor's memory back;
  // Eigen's solvers are noncopyable, so there is no way to reset one in place.
  using selinv_llt_t = Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int>>;
  std::unique_ptr<selinv_llt_t> selinv_llt;
  Eigen::SparseMatrix<double, 0, int> M_sym_;
  bool M_sym_ready{false};
  bool selinv_llt_ready{false};
  int selinv_viable_{-1}; // -1 undecided, 0 fill too high, 1 usable
  // Set once the caller has ruled the selected inverse out for this fit.
  // Everything the selinv path owns -- its private factorization, the
  // symmetric copy of the matrix it factorizes, and the selected inverse
  // itself -- is then dead weight, and for a 2-d mesh the private factor
  // alone is the same order as the solver's own. See disable_selinv().
  bool selinv_off_{false};
  void ensure_U(unsigned int seed);
  void ensure_QU(unsigned int seed);
  double reduce_probes(const Eigen::MatrixXd &MQU);
  // For non-symmetric mode we keep the last K to build normal equations and to
  // apply K^T on RHS when required
  Eigen::SparseMatrix<double, 0, int> K_last;

public:
  sparse_llt_solver() = default;
  sparse_llt_solver(int stype, int nin, int Ntrace, bool symmetric)
      : solver_type(ngme_fork_safe_stype(stype)), n(nin), N_iter(Ntrace),
        isSymmetric(symmetric) {}

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
    // Accelerate is not fork-safe; in a forked child this returns the CHOLMOD
    // equivalent instead. See ngme_fork_safe_stype() above.
    solver_type = ngme_fork_safe_stype(stype);
    isSymmetric = symmetric;
    // U / QU are the n x N_iter Hutchinson probe blocks. ensure_U() sizes them
    // on the first trace() / trace2() call, so a solver that is never asked for
    // a trace -- or a fit with Rao-Blackwellisation off -- never pays for them.
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
    {
      const int new_n = isSymmetric ? M.rows() : M.cols();
      if (new_n != n)
        U_computed = false; // probes are sized for the old dimension
      n = new_n;
    }
  }

  void compute(const Eigen::SparseMatrix<double, 0, int> &M) {
    if (!use_lu && isSymmetric && !selinv_off_) {
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
      U_computed = false; // probes are sized for the old dimension
    n = new_n;
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
  // tr(Q^-1 A^T diag(d) B) without ever forming A^T diag(d) B. The estimator
  // only ever needs that product applied to the probe block, and
  // (A^T diag(d) B) QU = A^T (d .* (B QU)) is three sparse-times-dense products
  // on an n x N_iter block instead of a sparse-sparse-sparse product whose
  // result is as dense as Q itself.
  double trace_factored(const Eigen::SparseMatrix<double, 0, int> &A,
                        const Eigen::VectorXd &d,
                        const Eigen::SparseMatrix<double, 0, int> &B,
                        unsigned int seed = 0);
  // Variance of the individual probe values in the last trace call (0 if exact).
  double last_probe_var() const { return last_probe_var_; }
  // Simplicial factorizations expose matrixL(); the CHOLMOD ones do not, and
  // the LU path is not a Cholesky at all.
  bool selinv_supported() const {
    return !use_lu && isSymmetric && !selinv_off_;
  }
  // Release the selected-inverse machinery for good. Called once the fill
  // test has come out against it, so neither the private factor nor the
  // per-compute() copy of the matrix is carried for the rest of the run.
  void disable_selinv() {
    selinv_off_ = true;
    selinv_llt_ready = false;
    S_sel_ready = false;
    M_sym_ready = false;
    selinv_llt.reset();
    S_sel = Eigen::SparseMatrix<double, 0, int>();
    M_sym_ = Eigen::SparseMatrix<double, 0, int>();
  }
  // nnz(L)/n. The cost of the selected inverse scales with the fill of the
  // factor, so this is what decides whether it beats probing.
  bool ensure_selinv_factor();
  // nnz(L)/n, which requires the factor and so, on a backend that does not
  // expose one, a private factorization of its own.
  double fill_ratio();
  // A lower bound on what fill_ratio() would return, read straight off the
  // matrix. The pattern of the Cholesky factor always contains the lower
  // triangle of the matrix it factorizes, so nnz(L) >= nnz(tril(M)). When even
  // that exceeds the threshold -- which it does for any 2-d spatial mesh -- the
  // selected inverse can be ruled out without factorizing anything. Returns 0
  // (i.e. rules nothing out) when the matrix is not held.
  double fill_lower_bound() const {
    if (!M_sym_ready || n <= 0)
      return 0.0;
    Eigen::Index cnt = 0;
    for (int c = 0; c < M_sym_.outerSize(); ++c)
      for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(M_sym_, c); it;
           ++it)
        if (it.row() >= c)
          ++cnt;
    return (double)cnt / (double)n;
  }
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
