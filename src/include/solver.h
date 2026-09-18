#ifndef __Solver__Solver__
#define __Solver__Solver__

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

#include "MatrixAlgebra.h"
#include "probing.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cholmod.h>
#include <atomic>
#include <cstdio>
#include <iostream>
#include <memory>
#include <utility>
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
// Fill-reducing ordering for the Accelerate backend, set once per fit from
// control_opt(solver_order=): -1 leaves the library default in place, -2 picks
// by fill, and a non-negative value is a SparseOrder_t. See analyze().
//
// Atomic only because it is written on the R thread at construction and read
// inside the OpenMP regions over chains and replicates; it never changes during
// a fit, so the relaxed ordering is all that is needed.
inline std::atomic<int> &ngme_accel_order_ref() {
  static std::atomic<int> v{-1};
  return v;
}
inline int ngme_accel_order() {
  return ngme_accel_order_ref().load(std::memory_order_relaxed);
}
inline void ngme_set_accel_order(int v) {
  ngme_accel_order_ref().store(v, std::memory_order_relaxed);
}
// ---------------------------------------------------------------------------
// Apple's Accelerate sparse solvers are NOT fork-safe.
//
// Eigen's AccelerateSupport wraps Apple's Sparse Solvers, which parallelise
// internally through libdispatch (GCD). GCD is documented as unusable in a
// process forked from one that has already initialised it: the child must
// exec() first. R's `parallel::mclapply()` forks and does not exec, so ANY use
// of an Accelerate factorization inside an mclapply worker aborts the child --
// and, because R's fork workers share the session's error handling, takes the
// whole R session down with it.
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
  // Ordering picked by fill under solver_order = "auto"; -1 until chosen.
  int accel_order_chosen_{-1};
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
  // Probe columns actually solved for. Equal to the budget below except under
  // structured probing, where it is quantized down to a whole number of
  // colourings, and on the exact path, where it is the dimension.
  int N_iter{10};
  // The budget in force: what the caller asked for, except that structured
  // probing may raise it to the smallest budget at which it can engage (see
  // probing_floor_). Kept apart from N_iter because the adaptive controller
  // steers this number and reads it back: a quantization that lands short must
  // not be read back as the new budget, or the budget could never climb past
  // the colour count it was rounded down to.
  int N_budget_{10};
  // What the caller last asked for, before any probing raise. Kept so that a
  // caller which re-asserts its own budget every iteration -- Operator does --
  // does not read the raised value back as a disagreement and undo the raise.
  int N_requested_{10};
  bool isSymmetric{true};
  Eigen::MatrixXd U, QU;
  bool QU_computed{false};
  // Hutchinson probe vectors.
  bool U_computed{false};
  unsigned int U_seed{0};
  // Set when U is the scaled identity rather than random probes, i.e. when the
  // probe budget is at least the dimension (n <= N_iter) and the trace is
  // therefore computed exactly.
  bool exact_trace{false};
  // Spread of the Hutchinson estimate from the most recent trace()/trace2()
  // call, scaled so that the variance the estimator contributes to the gradient
  // is last_probe_var_ / N_budget_ under EVERY probe scheme. This is what lets
  // the probe count be tuned against the Gibbs noise instead of being fixed a
  // priori. Zero when the trace was taken exactly, and -1 when probing is on
  // with a single replicate, where there is nothing to measure a spread from.
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
  // Set once the caller has ruled the selected inverse out for this fit.
  // Everything the selinv path owns -- its private factorization, the
  // symmetric copy of the matrix it factorizes, and the selected inverse
  // itself -- is then dead weight, and for a 2-d mesh the private factor
  // alone is the same order as the solver's own. See disable_selinv().
  bool selinv_off_{false};
  // --- structured probing (see include/probing.h) --------------------------
  // The index graph the probes are structured against, held only while probing
  // is a live possibility.
  ngme_probing::Adjacency probe_adj_;
  // Colourings already computed, by distance. Kept because the budget moves
  // during a fit and the distance is re-selected against it; each colouring is
  // then paid for once rather than at every budget change. A distance whose
  // count exceeded the cap is recorded as an empty colouring, and since the
  // count only grows with distance, no larger distance is tried after one.
  std::vector<std::vector<int>> colour_cache_;
  std::vector<int> colour_count_cache_;
  // The colouring in use, its colour count, the distance it came from, and how
  // many independent sign draws are averaged over it. While probing is on,
  // N_iter == n_colours_ * probe_reps_ <= N_budget_.
  //
  // The colouring is identified by its DISTANCE and looked up in the cache at
  // the point of use. Holding a pointer into colour_cache_ instead would not
  // survive the search in configure_probing(): testing the next distance grows
  // that cache, which moves everything already in it.
  int n_colours_{0};
  int probe_dist_{0};
  int probe_reps_{1};
  bool probing_{false};
  // Set by set_probing_source(); probing is then configured lazily, so a solver
  // that is never asked for a trace never pays for a colouring.
  bool probing_requested_{false};
  int probe_max_dist_{4};
  int probe_max_colours_{200};
  // Replicates required before probing is worth selecting. Two whenever the
  // probe variance has to be reported, because only whole replicates are
  // i.i.d. draws and a single one gives nothing to estimate a spread from.
  int probe_min_reps_{1};
  // Largest budget probing may raise itself to in order to engage at all.
  // Zero leaves the caller's budget alone, which is the historical behaviour
  // and what `trace_probing_raise_budget = 1` asks for.
  int probe_raise_cap_{0};
  bool probe_config_stale_{true};

  void ensure_U(unsigned int seed);
  void ensure_QU(unsigned int seed);
  double reduce_probes(const Eigen::MatrixXd &MQU);
  // Pick the colouring distance that fits the current budget, largest first,
  // and set N_iter from it. Cheap after the first call: the colourings are
  // cached and the budget rarely crosses a boundary.
  void configure_probing();
  // Ensures the colouring for distance d is cached; returns its colour count,
  // or 0 if it ran past the cap. The colouring itself is read out of the cache
  // by index, never held across another call.
  int colouring_for(int d);
  // The smallest budget at which probing engages, if the caller has allowed the
  // budget to be raised that far; 0 otherwise.
  int probing_floor_();
  // Recompute N_budget_ from the requested budget and the probing floor, and
  // invalidate whatever the old budget sized.
  void refresh_budget_();
  // For non-symmetric mode we keep the last K to build normal equations and to
  // apply K^T on RHS when required
  Eigen::SparseMatrix<double, 0, int> K_last;

public:
  sparse_llt_solver() = default;
  sparse_llt_solver(int stype, int nin, int Ntrace, bool symmetric)
      : solver_type(ngme_fork_safe_stype(stype)), n(nin), N_iter(Ntrace),
        N_budget_(Ntrace), N_requested_(Ntrace), isSymmetric(symmetric) {}

  // nonsym_mode: 0 = LU of K (default), 1 = Cholesky of the normal equations
  // K^T K (the historical path). Chosen by control_opt(nonsym_solver = ).
  inline void init(int nin, int Ntrace, bool symmetric, int stype,
                   int nonsym_mode = 0) {
    use_lu = !symmetric && nonsym_mode == 0;
    lu_ok = false;
    n = nin;
    N_iter = Ntrace;
    N_budget_ = Ntrace;
    N_requested_ = Ntrace;
    // Accelerate is not fork-safe; in a forked child this returns the CHOLMOD
    // equivalent instead. See ngme_fork_safe_stype() above.
    solver_type = ngme_fork_safe_stype(stype);
    isSymmetric = symmetric;
    // U / QU are the n x N_iter Hutchinson probe blocks. ensure_U() sizes them
    // on the first trace() / trace2() call, so a solver that is never asked for
    // a trace -- or a fit with Rao-Blackwellisation off -- never pays for them.
    QU_computed = false;
    S_sel_ready = false;
    U_computed = false;
    exact_trace = false;
    disable_probing();
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
    case 4: {
      // Accelerate's fill-reducing ordering, which the package never set: every
      // factorization used the library default (AMD). solver_order selects one
      // instead -- 2 = AMD, 3 = Metis nested dissection. Must be set BEFORE
      // analyzePattern.
      //
      // -2 chooses by FILL: the symbolic phase reports the size the numeric
      // factor will occupy, so the candidates are compared on that rather than
      // by timing a factorization.
      const Eigen::SparseMatrix<double, 0, int> Ana =
          isSymmetric ? M : (Eigen::SparseMatrix<double, 0, int>)(M.transpose() * M);
      if (ngme_accel_order() == -2 && accel_order_chosen_ < 0) {
        // Deliberately NOT including 1 (SparseOrderUser with a null
        // permutation, i.e. no permutation at all). It gives the SMALLEST
        // FACTOR of any candidate on a banded operator and is markedly SLOWER to
        // factorize regardless, because the factorization is supernodal and
        // the natural ordering of a band leaves supernodes too thin to reach
        // the blocked kernels. Fill ranks fill-reducing orderings against each
        // other reliably; it does not rank them against not ordering at all.
        static const int cand[] = {0, 2, 3};
        std::size_t best = 0;
        for (int c : cand) {
          R_accelerate.setOrder((SparseOrder_t)c);
          R_accelerate.analyzePattern(Ana);
          if (R_accelerate.info() != Eigen::Success)
            continue;
          const std::size_t fs = R_accelerate.factorSize();
          if (fs > 0 && (best == 0 || fs < best)) {
            best = fs;
            accel_order_chosen_ = c;
          }
        }
        if (accel_order_chosen_ < 0)
          accel_order_chosen_ = 0;
      }
      const int ord =
          (ngme_accel_order() == -2) ? accel_order_chosen_ : ngme_accel_order();
      if (ord >= 0)
        R_accelerate.setOrder((SparseOrder_t)ord);
      R_accelerate.analyzePattern(Ana);
      break;
    }
    case 5:
      if (ngme_accel_order() >= 0)
        R_accelerate_ldlt.setOrder((SparseOrder_t)ngme_accel_order());
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

  // The draw and the Rao-Blackwellised conditional mean share a factorization
  // and differ only by the perturbation of the right-hand side. A triangular
  // solve is bound by reading the factor, so the two together traverse it once
  // instead of twice. The arithmetic is that of the two calls it replaces: the
  // same solve() path, batched over columns.
  inline void rMVN_mean(const SparseMatrix<double, 0, int> &G,
                        const SparseMatrix<double, 0, int> &H,
                        Eigen::VectorXd &mu, Eigen::VectorXd &z1,
                        Eigen::VectorXd &z2, Eigen::VectorXd &draw,
                        Eigen::VectorXd &mean) {
    Eigen::MatrixXd rhs(mu.size(), 2);
    rhs.col(0) = G.transpose() * z1 + H.transpose() * z2 + mu;
    rhs.col(1) = mu;
    Eigen::MatrixXd sol = solve(rhs);
    draw = sol.col(0);
    mean = sol.col(1);
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
  // trace2 split in two. S = K^-1 A K^-1 U depends only on A, but the H_K block
  // wants tr(K^-1 dK_k K^-1 dK_j) over PAIRS, so doing the solve inside
  // trace2() repeats it n_theta_K(n_theta_K+1)/2 times where n_theta_K would
  // do. Identical arithmetic, with the invariant hoisted.
  Eigen::MatrixXd trace2_lhs(const Eigen::SparseMatrix<double, 0, int> &A,
                             unsigned int seed = 0);
  double trace2_reduce(const Eigen::SparseMatrix<double, 0, int> &B,
                       const Eigen::MatrixXd &S) const;
  double trace(const Eigen::SparseMatrix<double, 0, int> &,
               unsigned int seed = 0);
  // tr(Q^-1 A^T diag(d) B) without ever forming A^T diag(d) B. The estimator
  // only ever needs that product applied to the probe block, and
  // (A^T diag(d) B) QU = A^T (d .* (B QU)) is three sparse-times-dense products
  // on an n x N_iter block instead of a sparse-sparse-sparse product whose
  // result is as dense as Q itself.
  // trace_factored split so the B QU product can be shared. Both RB trace loops
  // pass the SAME B (= K), so recomputing B*QU inside every call repeats a
  // sparse-times-dense product once per parameter for no reason. d still varies
  // per parameter, but scaling the shared block by it is elementwise.
  Eigen::MatrixXd trace_factored_rhs(const Eigen::SparseMatrix<double, 0, int> &B,
                                     unsigned int seed = 0);
  double trace_factored_with(const Eigen::SparseMatrix<double, 0, int> &A,
                             const Eigen::VectorXd &d,
                             const Eigen::MatrixXd &BQU);
  double trace_factored(const Eigen::SparseMatrix<double, 0, int> &A,
                        const Eigen::VectorXd &d,
                        const Eigen::SparseMatrix<double, 0, int> &B,
                        unsigned int seed = 0);
  // Spread of the last trace estimate, scaled so that the estimator's variance
  // is this over get_N_iter(). Zero if the trace was exact, negative if it
  // could not be measured (see last_probe_var_).
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
  // Stored entries of the factor and of the selected inverse. These are the
  // deterministic counterparts of "how long does each route take": a probe is
  // a triangular solve, costing nnz(L), and the selected inverse is built and
  // then read, costing its own nnz. Using sizes rather than elapsed time keeps
  // the choice between them a function of the matrices alone, so a fit does not
  // depend on how busy the machine was.
  long long factor_nnz();
  // Flop count of the Takahashi recursion, read off the factor's pattern
  // without forming anything. Column j scatters its subdiagonal rows and then
  // walks the column of L belonging to each of them, so the work is
  //   sum_j sum_{r in L[:,j], r > j} nnz(L[:,r]),
  // which is the recursion's own operation count rather than a proxy for it.
  // nnz(S_sel) is not that count and cannot stand in for it: the selected
  // inverse shares L's pattern exactly, so its size says what the result costs
  // to store and nothing about what it costs to form. Sizing the build by it
  // makes the exact route look free on every matrix.
  long long selinv_build_flops();
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
  // The probe budget, which is what the adaptive controller steers and what
  // last_probe_var() is expressed against. Identical to the number of probe
  // columns unless structured probing quantized it down; see N_budget_.
  int get_N_iter() const { return N_budget_; }
  // What the caller last asked for, before any probing raise. A caller that
  // re-asserts its own budget each iteration must compare against THIS, not
  // against get_N_iter(), or it will read the raise back as a disagreement and
  // undo it on every pass.
  int get_requested_N_iter() const { return N_requested_; }
  // Probe columns actually solved for.
  int get_n_probes() const { return N_iter; }
  bool is_exact_trace() const { return exact_trace; }
  // Change the probe budget at run time. Invalidates the cached probes so the
  // next trace call regenerates them (and switches to the exact path if the
  // budget now reaches the dimension).
  void set_N_iter(int Ntrace) {
    if (Ntrace < 1 || Ntrace == N_requested_)
      return;
    N_requested_ = Ntrace;
    refresh_budget_();
  }

  // Structure the Hutchinson probes against the graph of A instead of drawing
  // them densely. Costs nothing here: the adjacency is built, but the colouring
  // itself waits until a trace is actually asked for, and is then chosen
  // against the budget in force at that moment. Passing max_dist < 1 -- or a
  // matrix whose graph cannot be coloured inside max_colours -- leaves the
  // solver on plain Rademacher probes.
  //
  // min_reps is the number of independent sign draws the caller needs over one
  // colouring. Two are required to report a probe variance, since only whole
  // replicates are exchangeable; the distance is then chosen so that min_reps
  // of them still fit the budget, which is what keeps the probe count from
  // rising above what was already being spent.
  //
  // raise_cap is the largest budget probing may raise ITSELF to when the
  // caller's budget is too small for any colouring. Zero leaves the budget
  // alone, and probing then simply stays off at such a budget.
  void set_probing_source(const Eigen::SparseMatrix<double, 0, int> &A,
                          int max_dist, int max_colours, int min_reps,
                          int raise_cap = 0) {
    disable_probing();
    if (max_dist < 1 || A.rows() != A.cols() || A.rows() <= 0)
      return;
    probe_adj_ = ngme_probing::build_adjacency(A);
    if (probe_adj_.empty())
      return;
    probe_max_dist_ = max_dist;
    probe_max_colours_ = std::max(1, max_colours);
    probe_min_reps_ = std::max(1, min_reps);
    probe_raise_cap_ = std::max(0, raise_cap);
    probing_requested_ = true;
    refresh_budget_();
  }
  void disable_probing() {
    const bool was_raised = N_budget_ != N_requested_;
    probing_requested_ = false;
    probing_ = false;
    probe_config_stale_ = true;
    n_colours_ = 0;
    probe_dist_ = 0;
    probe_reps_ = 1;
    probe_raise_cap_ = 0;
    colour_cache_.clear();
    colour_count_cache_.clear();
    probe_adj_.clear();
    // A budget raised to reach a colouring goes back to what was asked for.
    if (was_raised)
      refresh_budget_();
    else
      N_budget_ = N_requested_;
  }
  // Smallest budget at which structured probing would engage, or 0 if it never
  // would. The budget controller uses this as a floor once probing is already
  // on; see the note at its call site for why only then.
  int probe_min_budget() {
    if (!probing_requested_)
      return 0;
    const int p = colouring_for(1);
    return p > 0 ? p * probe_min_reps_ : 0;
  }
  bool probing_active() const { return probing_; }
  int probe_colours() const { return n_colours_; }
  int probe_dist() const { return probe_dist_; }
  int probe_reps() const { return probe_reps_; }

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
