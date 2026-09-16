#include "operator.h"
#include <sstream>
#include "include/phase_timing.h"
#include <cstdio>
#include <cstdlib>
#include "include/factor_counters.h"
#include <algorithm>
#include <cmath>

// for initialize Latent models
Operator::~Operator() = default;

void Operator::record_K_pattern() {
  const int outer = K.outerSize() + 1;
  const int nnz = static_cast<int>(K.nonZeros());
  K_pat_outer.assign(K.outerIndexPtr(), K.outerIndexPtr() + outer);
  K_pat_inner.assign(K.innerIndexPtr(), K.innerIndexPtr() + nnz);
  analyzed_cholK = true;
}

bool Operator::K_pattern_changed() const {
  if (!analyzed_cholK)
    return true;
  if (K_pat_outer.size() != static_cast<size_t>(K.outerSize()) + 1)
    return true;
  if (K_pat_inner.size() != static_cast<size_t>(K.nonZeros()))
    return true;
  if (!std::equal(K_pat_outer.begin(), K_pat_outer.end(), K.outerIndexPtr()))
    return true;
  return !std::equal(K_pat_inner.begin(), K_pat_inner.end(),
                     K.innerIndexPtr());
}

// Exact traces for a triangular operator, in O(n) and with no factorization.
//
// If K is triangular then so is K^{-1}, and for lower-triangular K and dK the
// only term surviving on the diagonal of the product is k = i:
//     (K^{-1} dK)_ii = sum_k (K^{-1})_ik (dK)_ki = (dK)_ii / K_ii,
// since (K^{-1})_ik = 0 for k > i and (dK)_ki = 0 for k < i (and symmetrically
// for upper-triangular). Hence
//     tr(K^{-1} dK)              = sum_i (dK)_ii / K_ii
//     tr(K^{-1} dK_k K^{-1} dK_j) = sum_i (dK_k)_ii (dK_j)_ii / K_ii^2
bool Operator::try_triangular_traces(const UpdateOptions &opts, bool want_trace,
                                     bool want_HK) {
  const int n = static_cast<int>(K.rows());
  if (n == 0 || K.rows() != K.cols())
    return false;

  // Orientations a matrix is compatible with: bit 0 = no entries above the
  // diagonal, bit 1 = none below. A diagonal matrix sets both.
  auto orient = [](const SparseMatrix<double> &M, int n) {
    if (M.rows() != n || M.cols() != n)
      return 0;
    int o = 3;
    for (int c = 0; c < M.outerSize() && o; ++c) {
      for (SparseMatrix<double>::InnerIterator it(M, c); it; ++it) {
        if (it.value() == 0.0)
          continue;
        if (it.row() < it.col())
          o &= ~1;
        if (it.row() > it.col())
          o &= ~2;
        if (!o)
          break;
      }
    }
    return o;
  };

  // The cancellation that makes the diagonal formula exact needs K and every
  // derivative used to share one orientation: (K^{-1})_ik = 0 for k > i kills
  // the sum only when (dK)_ki = 0 for k < i as well. Checking K alone would be
  // wrong for an operator whose K happens to be diagonal at the current theta
  // while dK is not triangular.
  int o = orient(K, n);
  if (!o)
    return false;
  for (int j = 0; j < n_theta_K && o; ++j) {
    if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j])
      continue;
    if ((int)dK.size() <= j)
      return false;
    o &= orient(dK[j], n);
  }
  if (want_HK && !d2K.empty()) {
    for (int j = 0; j < n_theta_K && o; ++j)
      for (int k = j; k < n_theta_K && o; ++k)
        if (d2K[j][k].rows() == n)
          o &= orient(d2K[j][k], n);
  }
  if (!o)
    return false;

  VectorXd d(n);
  for (int i = 0; i < n; ++i) {
    const double v = K.coeff(i, i);
    if (!(std::abs(v) > 0.0))
      return false; // singular; let the generic path deal with it
    d(i) = v;
  }

  auto diag_ratio = [&](const SparseMatrix<double> &M) {
    double t = 0.0;
    for (int i = 0; i < n; ++i)
      t += M.coeff(i, i) / d(i);
    return t;
  };

  if (want_trace) {
    if (trace_vals.size() != n_theta_K)
      trace_vals = VectorXd::Zero(n_theta_K);
    for (int j = 0; j < n_theta_K; ++j) {
      if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j]) {
        trace_vals(j) = 0.0;
        continue;
      }
      trace_vals(j) = diag_ratio(dK[j]);
    }
  }

  if (want_HK) {
    if (HK_trace.rows() != n_theta_K || HK_trace.cols() != n_theta_K)
      HK_trace = MatrixXd::Zero(n_theta_K, n_theta_K);
    else
      HK_trace.setZero();
    MatrixXd r(n, n_theta_K); // (dK_j)_ii / K_ii
    for (int j = 0; j < n_theta_K; ++j)
      for (int i = 0; i < n; ++i)
        r(i, j) = dK[j].coeff(i, i) / d(i);
    for (int j = 0; j < n_theta_K; ++j) {
      for (int k = j; k < n_theta_K; ++k) {
        double t1 = 0.0;
        for (int i = 0; i < n; ++i)
          t1 -= r(i, k) * r(i, j);
        double t2 = 0.0;
        if (!d2K.empty() && d2K[j][k].rows() == n)
          t2 = diag_ratio(d2K[j][k]);
        HK_trace(j, k) = t1 + t2;
        if (j != k)
          HK_trace(k, j) = HK_trace(j, k);
      }
    }
  }
  return true;
}

namespace {
// Same sparsity pattern? Two builds of one operator at perturbed theta differ
// only in their values, so this is true on every iteration after the first.
template <class M> bool same_pattern(const M &a, const M &b) {
  if (a.rows() != b.rows() || a.cols() != b.cols() ||
      a.nonZeros() != b.nonZeros() || !a.isCompressed() || !b.isCompressed())
    return false;
  return std::equal(a.outerIndexPtr(), a.outerIndexPtr() + a.outerSize() + 1,
                    b.outerIndexPtr()) &&
         std::equal(a.innerIndexPtr(), a.innerIndexPtr() + a.nonZeros(),
                    b.innerIndexPtr());
}

// out = (x - y) * s, reusing out's storage when all three already share one
// pattern. Eigen's sparse subtract merges patterns and allocates a fresh matrix
// every call; here the operands are the same operator at perturbed parameters,
// so the pattern never moves and a single pass over the values suffices. This
// sits inside the differencing loop, which is the dominant cost for an operator
// without analytic derivatives.
template <class M>
void diff_into(M &out, const M &x, const M &y, double s) {
  if (same_pattern(x, y) && same_pattern(out, x)) {
    const double *xv = x.valuePtr(), *yv = y.valuePtr();
    double *ov = out.valuePtr();
    const Eigen::Index nnz = x.nonZeros();
    for (Eigen::Index i = 0; i < nnz; ++i)
      ov[i] = (xv[i] - yv[i]) * s;
    return;
  }
  out = (x - y) * s; // establishes the pattern; the fast path takes over next
  out.makeCompressed();
}
} // namespace

// Unified updater: update K, Z, (optionally) dK, dZ, factorization, and traces
void Operator::update_all(const VectorXd &theta, const UpdateOptions &opts) {
  // Visible to update_dKdZ / update_d2Kd2Z, whose signatures carry no options.
  // Cleared on the way out so a stale pointer cannot be read later.
  cur_opts_ = &opts;
  struct OptsGuard {
    const UpdateOptions **slot;
    ~OptsGuard() { *slot = nullptr; }
  } _og{&cur_opts_};

  // Every probe drawn anywhere below this point is solved against K, not
  // against the block precision -- this class holds no other solver. Marking
  // the whole call rather than just the triangular-trace helper is the point:
  // the general trace path further down (and the H_K block) probe outside that
  // helper, so with the marker set only there those solves were attributed to
  // the block-precision counter and k_probe_solves stayed at zero however hard
  // the operator was probing.
  ngme_counters::probe_role_scope _role(ngme_counters::probe_role::op);

  // 1) Build K and Z once at base theta
  { ngme_timing::Scope _s(ngme_timing::op_build_us()); build_KZ(theta); }

  // Ensure storage ready
  if ((int)dK.size() != n_theta_K)
    dK.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));
  if ((int)dZ.size() != n_theta_K)
    dZ.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));
  // d2K / d2Z are n_theta_K^2 matrices of size n x n. Even empty, each one
  // carries an outer-index array of n + 1 ints, so allocating them for an
  // operator that is never asked for second derivatives costs
  // 4 * n_theta_K^2 * n bytes for nothing. Allocate on demand instead.
  const bool want_d2 = opts.compute_d2K || opts.compute_d2Z;
  if (want_d2 && (int)d2K.size() != n_theta_K)
    d2K.assign(n_theta_K,
               std::vector<SparseMatrix<double>>(
                   n_theta_K, SparseMatrix<double>(h.size(), h.size())));
  if (want_d2 && (int)d2Z.size() != n_theta_K)
    d2Z.assign(n_theta_K, std::vector<SparseMatrix<double, 0, int>>(
                              n_theta_K, SparseMatrix<double, 0, int>(
                                             h.size(), h.size())));

  // 2) Derivatives: analytic or numeric
  ngme_timing::Scope *_dk = new ngme_timing::Scope(ngme_timing::op_dK_us());
  bool have_analytic = false;
  last_eps_dK_ = opts.eps_dK;
  if ((opts.compute_dK || opts.compute_dZ) && opts.prefer_analytic_dK) {
    have_analytic = update_dKdZ(theta);
  }
  if ((opts.compute_dK || opts.compute_dZ) && !have_analytic) {
    // Numeric differencing around theta
    // Backup base K/Z
    SparseMatrix<double> K_base = K;
    SparseMatrix<double> Z_base = Z;
    auto is_fixed = [&](int j) {
      return (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j]);
    };
    for (int j = 0; j < n_theta_K; ++j) {
      if (is_fixed(j)) {
        dK[j].setZero();
        dZ[j].setZero();
        continue;
      }
      double eps = opts.eps_dK;
      if (opts.diff_dK_mode == DiffMode::Forward) {
        VectorXd th_f = theta;
        th_f(j) += eps;
        build_KZ(th_f);
        if (opts.compute_dK)
          diff_into(dK[j], K, K_base, 1.0 / eps);
        if (opts.compute_dZ)
          diff_into(dZ[j], Z, Z_base, 1.0 / eps);
      } else {
        VectorXd th_p = theta;
        th_p(j) += eps;
        build_KZ(th_p);
        SparseMatrix<double> Kp = K;
        SparseMatrix<double, 0, int> Zp = Z;
        VectorXd th_m = theta;
        th_m(j) -= eps;
        build_KZ(th_m);
        SparseMatrix<double> Km = K;
        SparseMatrix<double, 0, int> Zm = Z;
        if (opts.compute_dK)
          diff_into(dK[j], Kp, Km, 1.0 / (2.0 * eps));
        if (opts.compute_dZ)
          diff_into(dZ[j], Zp, Zm, 1.0 / (2.0 * eps));
      }
      // Restore base
      K = K_base;
      Z = Z_base;
    }
  }
  // 2b) Second derivatives (numeric unless analytic provided)
  bool have_analytic2 = false;
  if ((opts.compute_d2K || opts.compute_d2Z) && opts.prefer_analytic_d2K) {
    have_analytic2 = update_d2Kd2Z(theta);
  }
  if ((opts.compute_d2K || opts.compute_d2Z) && !have_analytic2) {
    // Mixed partials central difference for j!=k, and pure second for j==k
    SparseMatrix<double> K0 = K;
    SparseMatrix<double, 0, int> Z0 = Z;
    auto step_for = [&](double th) {
      double s = std::max(1.0, std::abs(th));
      return opts.eps_dK * s;
    };
    for (int j = 0; j < n_theta_K; ++j) {
      double hj = step_for(theta(j));
      // Pure second K_{jj}
      if (opts.compute_d2K) {
        VectorXd thp = theta;
        thp(j) += hj;
        build_KZ(thp);
        SparseMatrix<double> Kp = K;
        VectorXd thm = theta;
        thm(j) -= hj;
        build_KZ(thm);
        SparseMatrix<double> Km = K;
        d2K[j][j] = (Kp - 2.0 * K0 + Km) * (1.0 / (hj * hj));
      }
      if (opts.compute_d2Z) {
        VectorXd thp = theta;
        thp(j) += hj;
        build_KZ(thp);
        auto Zp = Z;
        VectorXd thm = theta;
        thm(j) -= hj;
        build_KZ(thm);
        auto Zm = Z;
        d2Z[j][j] = (Zp - 2.0 * Z0 + Zm) * (1.0 / (hj * hj));
      }
      // Mixed K_{jk}, j<k
      for (int k = j + 1; k < n_theta_K; ++k) {
        double hk = step_for(theta(k));
        VectorXd th_pp = theta;
        th_pp(j) += hj;
        th_pp(k) += hk;
        build_KZ(th_pp);
        auto Kpp = K;
        auto Zpp = Z;
        VectorXd th_pm = theta;
        th_pm(j) += hj;
        th_pm(k) -= hk;
        build_KZ(th_pm);
        auto Kpm = K;
        auto Zpm = Z;
        VectorXd th_mp = theta;
        th_mp(j) -= hj;
        th_mp(k) += hk;
        build_KZ(th_mp);
        auto Kmp = K;
        auto Zmp = Z;
        VectorXd th_mm = theta;
        th_mm(j) -= hj;
        th_mm(k) -= hk;
        build_KZ(th_mm);
        auto Kmm = K;
        auto Zmm = Z;
        if (opts.compute_d2K) {
          d2K[j][k] = (Kpp - Kpm - Kmp + Kmm) * (1.0 / (4.0 * hj * hk));
          d2K[k][j] = d2K[j][k];
        }
        if (opts.compute_d2Z) {
          auto Zmix = (Zpp - Zpm - Zmp + Zmm) * (1.0 / (4.0 * hj * hk));
          d2Z[j][k] = Zmix;
          d2Z[k][j] = Zmix;
        }
      }
      // Restore base for next j
      K = K0;
      Z = Z0;
    }
  }

  // NGME_DIFF_CHECK=1: whenever a closed form supplied dK or d2K, recompute the
  // same quantity by differencing and report the relative disagreement at two
  // step sizes. Here the NUMERIC value is the approximation, so a correct
  // closed form makes the error SHRINK with the step, roughly linearly for
  // the forward dK, quadratically for the central d2K. A wrong closed form
  // leaves the error flat, which is what makes this a test rather than a
  // tolerance check.
  if ((have_analytic || have_analytic2) && std::getenv("NGME_DIFF_CHECK")) {
    const SparseMatrix<double> K_save = K;
    const SparseMatrix<double, 0, int> Z_save = Z;
    auto relerr = [](const SparseMatrix<double> &x,
                     const SparseMatrix<double> &y) {
      const double d = (SparseMatrix<double>(x - y)).norm();
      const double n = std::max(1e-300, y.norm());
      return d / n;
    };
    // Step sizes to sweep. The d2K reference divides by 4 h^2, so at the
    // shipped eps=1e-4 that is 4e-8 and rounding noise in K swamps the
    // truncation error -- the check has to run where truncation dominates, and
    // then watch it fall like h^2. Override with NGME_DIFF_CHECK_EPS.
    std::vector<double> eps_list{1e-2, 3e-3, 1e-3};
    if (const char *el = std::getenv("NGME_DIFF_CHECK_EPS")) {
      eps_list.clear();
      std::string sl(el), tok;
      std::stringstream ss(sl);
      while (std::getline(ss, tok, ','))
        if (!tok.empty()) eps_list.push_back(std::atof(tok.c_str()));
    }
    for (size_t pass = 0; pass < eps_list.size(); ++pass) {
      const double e = eps_list[pass];
      double worst1 = 0.0, worst2 = 0.0;
      int w1 = -1, w2j = -1, w2k = -1;
      if (have_analytic)
        for (int j = 0; j < n_theta_K; ++j) {
          VectorXd th = theta; th(j) += e;
          build_KZ(th);
          SparseMatrix<double> ref = (K - K_save) * (1.0 / e);
          const double r = relerr(dK[j], ref);
          if (r > worst1) { worst1 = r; w1 = j; }
        }
      if (have_analytic2)
        for (int j = 0; j < n_theta_K; ++j)
          for (int k = j; k < n_theta_K; ++k) {
            const double hj = e * std::max(1.0, std::abs(theta(j)));
            const double hk = e * std::max(1.0, std::abs(theta(k)));
            SparseMatrix<double> ref;
            if (j == k) {
              VectorXd tp2 = theta; tp2(j) += hj; build_KZ(tp2);
              SparseMatrix<double> Kp = K;
              VectorXd tm2 = theta; tm2(j) -= hj; build_KZ(tm2);
              ref = (Kp - 2.0 * K_save + K) * (1.0 / (hj * hj));
            } else {
              VectorXd t1 = theta; t1(j) += hj; t1(k) += hk; build_KZ(t1);
              SparseMatrix<double> Kpp = K;
              VectorXd t2 = theta; t2(j) += hj; t2(k) -= hk; build_KZ(t2);
              SparseMatrix<double> Kpm = K;
              VectorXd t3 = theta; t3(j) -= hj; t3(k) += hk; build_KZ(t3);
              SparseMatrix<double> Kmp = K;
              VectorXd t4 = theta; t4(j) -= hj; t4(k) -= hk; build_KZ(t4);
              ref = (Kpp - Kpm - Kmp + K) * (1.0 / (4.0 * hj * hk));
            }
            const double r = relerr(d2K[j][k], ref);
            if (r > worst2) { worst2 = r; w2j = j; w2k = k; }
          }
      std::fprintf(stderr,
                   "[diff check] eps=%.3g  dK worst %.3e (theta_K[%d])  "
                   "d2K worst %.3e (%d,%d)\n",
                   e, worst1, w1, worst2, w2j, w2k);
    }
    K = K_save;
    Z = Z_save;
  }

  delete _dk; // closes the dK timer
  ngme_timing::Scope _tr(ngme_timing::op_trace_us());
  // 3) Factorization. cholK_solver is used for nothing but the traces below,
  // so it is only worth building when a trace is actually wanted and no
  // structural shortcut supplies it.
  trace_ready = false;
  HK_trace_ready = false;
  const bool want_trace =
      opts.compute_trace && !dK.empty() && n_theta_K > 0;
  const bool want_HK = opts.compute_HK_trace && n_theta_K > 0;
  if (!want_trace && !want_HK)
    return;
  if (try_triangular_traces(opts, want_trace, want_HK)) {
    trace_ready = want_trace;
    HK_trace_ready = want_HK;
    return;
  }
  if (compute_traces_structured(theta, opts)) {
    trace_ready = want_trace;
    HK_trace_ready = want_HK;
    return;
  }

  if (!llt_inited) {
    // If opts.solver_type not set by caller, stays default 0
    cholK_solver.init(K.rows(), opts.n_trace_iter, symmetric, opts.solver_type,
                      opts.nonsym_solver);
    llt_inited = true;
  } else if (opts.n_trace_iter > 0 &&
             cholK_solver.get_requested_N_iter() != opts.n_trace_iter) {
    // The budget is initialised once but may be driven at run time. Picking it
    // up here, before any trace is taken, keeps the operator-side probe count
    // in step with the caller instead of frozen at its construction value --
    // and these probes feed the gradient of the operator parameters directly.
    cholK_solver.set_N_iter(opts.n_trace_iter);
  }
  // Storage must be packed before the pattern can be compared (and before the
  // solvers see it).
  if (!K.isCompressed())
    K.makeCompressed();
  // Re-run the symbolic phase exactly when the sparsity pattern moved.
  const bool K_pattern_moved = K_pattern_changed();
  if (K_pattern_moved || ngme_counters::cache_disabled()) {
    ngme_counters::bump(ngme_counters::K_analyzes);
    { ngme_timing::Scope _s(ngme_timing::k_symbolic_us()); cholK_solver.analyze(K); }
    record_K_pattern();
  }
  // Re-source the colouring when the pattern REALLY moved, which is not the
  // same condition as re-running the symbolic phase. The cache-disabled
  // diagnostic re-runs that phase on an UNCHANGED pattern, to show the cache
  // is not hiding anything, and so must not change any number. Re-sourcing
  // here would change them: the colouring is rebuilt, and the limit on how
  // often it may be rebuilt is then reached in that mode alone, leaving the
  // two runs on different probe schemes. Only the graph is built here; the
  // colouring itself waits until a trace is asked for.
  if (K_pattern_moved)
    setup_K_probing(opts);
  { ngme_timing::Scope _s(ngme_timing::k_numeric_us()); cholK_solver.compute(K); }
  // Either route to K^{-1} for a non-symmetric operator -- the LU of K, or the
  // Cholesky of K^T K -- is factorizing something that is singular only if K
  // itself is degenerate at the current parameters, so a failure means the same
  // thing in both and the traces below would be meaningless either way. It is
  // also unsafe in the LU case: Eigen returns from factorize() before setting
  // up its L store and guards the solve path only with eigen_assert, which is
  // compiled out here, so continuing would read uninitialised memory. Report it
  // as an ordinary error.
  //
  // A symmetric operator keeps its historical behaviour: a failed Cholesky
  // there is memory-safe, and during SGD it is usually a transient excursion
  // the next iterate recovers from.
  if (!symmetric && !cholK_solver.factorization_success()) {
    std::string msg = "Factorization of the operator matrix K failed for model '" +
                      generic_type +
                      "': K is singular or numerically rank-deficient at the "
                      "current parameter values.";
    if (!cholK_solver.uses_lu())
      msg += " If K is merely ill-conditioned rather than singular, "
             "control_opt(nonsym_solver = \"lu\") factorizes K directly "
             "instead of K^T K and does not square its condition number.";
    throw std::runtime_error(msg);
  }

  // 4) Traces (only if dK available)
  if (want_trace) {
    if (trace_vals.size() != n_theta_K)
      trace_vals = VectorXd::Zero(n_theta_K);
    // tr(K^{-1} dK)
    for (int j = 0; j < n_theta_K; ++j) {
      if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j]) {
        trace_vals(j) = 0.0;
        continue;
      }
      // Exact trace from the selected inverse, or Hutchinson probes, decided
      // once per fit and latched.
      if (selinv_state_ < 0) {
        // The operator is decided by COUNTING the two routes against each
        // other, where BlockModel::qq_trace is decided by a fill threshold.
        // The two sides are not symmetric. Counting only settles the question
        // when it is not close, because the routes reach very different
        // fractions of peak throughput per operation -- a probe is a triangular
        // solve over a block of right-hand sides and vectorises better the
        // denser the factor gets, while the Takahashi recursion is a scalar
        // scatter/gather. On the block precision that spread is wider than the
        // margin between the models, so a count cannot arbitrate and the fill
        // threshold stays. Here it is not close: the operator's factor is far
        // sparser than the probes it replaces, its traces are taken once per
        // optimizer iteration rather than once per Gibbs pass, and for an
        // operator that factors into a tensor product the exact route
        // decomposes through the factors. Every model tested clears the
        // threshold by more than an order of magnitude.
        //
        // The free fill bound is kept ahead of the count purely as a rail: the
        // count needs the factor, and forming one on a hopeless matrix is the
        // cost the rail exists to avoid.
        const double max_fill = opts.selinv_max_fill;
        double fr = cholK_solver.fill_lower_bound();
        bool ruled_out =
            !cholK_solver.selinv_supported() || fr > max_fill;
        if (!ruled_out) {
          // Operation counts rather than timings: the choice must be a
          // function of the matrices, not of how busy the machine is, or the
          // same fit stops being reproducible -- and a reproducibility test
          // asserts exactly that. The exact
          // route is charged the Takahashi recursion that forms the selected
          // inverse, not the size of the result -- the two are unrelated, since
          // the selected inverse carries the factor's pattern whatever it cost
          // to build.
          const long long nnz_L = cholK_solver.factor_nnz();
          const long long build = cholK_solver.selinv_build_flops();
          const bool sel_ok = (nnz_L > 0 && build > 0);
          const double m = std::max(1.0, (double)n_theta_K);
          const long long dim = std::max<long long>(1, (long long)K.rows());
          const double reuse = (double)nnz_L / (double)dim;
          const double t_sel = (double)build + m * reuse;
          const double N = (double)std::max(1, cholK_solver.get_N_iter());
          // One probe block serves every parameter, so the solves are charged
          // once rather than once per parameter; a probe is a forward and a
          // back substitution, hence the 2.
          const double t_probe = 2.0 * N * (double)nnz_L + N * m * reuse;
          ruled_out = !(sel_ok && t_probe > 0.0 &&
                        t_sel <= opts.selinv_cost_ratio * t_probe);
          if (opts.debug)
            Rprintf("[selinv K] nnz_L=%lld build=%lld m=%g N=%g exact=%g "
                    "probes=%g\n",
                    nnz_L, build, m, N, t_sel, t_probe);
        }
        selinv_state_ = ruled_out ? 0 : 1;
        if (opts.debug)
          Rprintf("[selinv K] fill>=%g -> %s\n", fr,
                  selinv_state_ ? "exact selected inverse"
                                : "Hutchinson probes");
        if (selinv_state_ == 0)
          cholK_solver.disable_selinv(); // hand back the factor and S_sel
      }
      if (selinv_state_ == 1) {
        double v = 0.0;
        if (cholK_solver.selinv_trace(dK[j], v)) {
          trace_vals(j) = v;
          continue;
        }
        // dK reached outside the factor pattern; probes for the rest of the fit.
        selinv_state_ = 0;
        cholK_solver.disable_selinv();
      }
      trace_vals(j) = cholK_solver.trace(dK[j], opts.trace_seed);
    }
  }
  trace_ready = want_trace;

  // 5) H_K trace block: -tr(K^{-1}K_k K^{-1}K_j) + tr(K^{-1}K_{jk})
  if (want_HK) {
    if (HK_trace.rows() != n_theta_K || HK_trace.cols() != n_theta_K)
      HK_trace = MatrixXd::Zero(n_theta_K, n_theta_K);
    else
      HK_trace.setZero();
    // Loop over k OUTSIDE: S = K^-1 dK_k K^-1 U depends only on k, so it is
    // computed once per parameter instead of once per pair. That turns the
    // expensive solve from O(n_theta_K^2) into O(n_theta_K).
    for (int k = 0; k < n_theta_K; ++k) {
      Eigen::MatrixXd Sk;
      bool Sk_ready = false;
      for (int j = 0; j <= k; ++j) {
        double t1 = 0.0, t2 = 0.0;
        // -tr(K^{-1} K_k K^{-1} K_j)
        if (dK[k].rows() > 0 && dK[j].rows() > 0) {
          if (!Sk_ready) {
            Sk = cholK_solver.trace2_lhs(dK[k], opts.trace_seed);
            Sk_ready = true;
          }
          t1 = -cholK_solver.trace2_reduce(dK[j], Sk);
        }
        // + tr(K^{-1} K_{jk})
        if (!d2K.empty() && d2K[j][k].rows() > 0) {
          t2 = cholK_solver.trace(d2K[j][k], opts.trace_seed);
        }
        double val = t1 + t2;
        HK_trace(j, k) = val;
        if (j != k)
          HK_trace(k, j) = val;
      }
    }
    HK_trace_ready = true;
  }
}

// Point cholK_solver's probes at the graph of K, or take them off it.
//
// Both operator-side estimators ride the same probe block: tr(K^{-1} dK), where
// the colouring pays for itself, and the H_K pair tr(K^{-1}K_k K^{-1}K_j),
// where it does not -- two inverses roughly double the decay length, so the
// colouring distances that fit a probe budget barely reach it. Sharing costs
// nothing either way: H_K measures no worse on structured probes than on dense
// ones, and splitting the two apart would mean two probe blocks and two sets of
// solves where there is now one.
void Operator::setup_K_probing(const UpdateOptions &opts) {
  if (!opts.trace_probing) {
    cholK_solver.disable_probing();
    return;
  }
  // The colouring describes one sparsity pattern and is amortized over the
  // iterations that share it. An operator whose pattern keeps moving never
  // amortizes it, so probing is dropped rather than re-paid every iteration;
  // the same reasoning as BlockModel::setup_qq_probing().
  if (++k_probing_setups_ > 3) {
    cholK_solver.disable_probing();
    return;
  }
  const int budget = cholK_solver.get_requested_N_iter();
  const int max_colours = std::max(budget, opts.trace_probing_max_colours);
  // Two caps on raising the budget to reach a colouring, tighter wins: a
  // multiple of what the caller asked for, and the ceiling this fit may spend.
  const int raise_cap =
      opts.trace_probing_raise_budget > 1.0
          ? std::min((int)std::floor(opts.trace_probing_raise_budget * budget),
                     std::max(budget, opts.trace_probing_max_colours))
          : 0;
  // One sign draw is enough here: nothing reads the operator-side probe
  // variance, so there is no spread to measure and no reason to spend a second
  // replicate on one.
  cholK_solver.set_probing_source(K, opts.trace_probing_max_dist, max_colours,
                                  /*min_reps*/ 1, raise_cap);
}

std::shared_ptr<Operator>
OperatorFactory::create(const Rcpp::List &operator_in) {
  string model_type = Rcpp::as<string>(operator_in["model"]);
  string generic_type = Rcpp::as<string>(operator_in["generic_type"]);

  if (model_type == "generic") {
    return std::make_shared<Generic>(operator_in);
  } else if (model_type == "generic_ns") {
    return std::make_shared<generic_ns>(operator_in);
  } else if (generic_type == "generic") {
    return std::make_shared<Generic>(operator_in);
  } else if (generic_type == "generic_ns") {
    return std::make_shared<generic_ns>(operator_in);
  } else if (model_type == "tp") {
    return std::make_shared<Tensor_prod>(operator_in);
  } else if (model_type == "spacetime") {
    return std::make_shared<Spacetime>(operator_in);
  } else if (model_type == "matern") {
    // Unify: always use the new Matern operator.
    // It handles stationary/non-stationary, integer/fractional alpha, and
    // produces Z.
    return std::make_shared<Matern>(operator_in);
  } else if (model_type == "arma") {
    return std::make_shared<Arma>(operator_in);
  } else if (model_type == "ou") {
    return std::make_shared<OU>(operator_in);
  } else if (model_type == "re") {
    return std::make_shared<Randeff>(operator_in);
  } else if (model_type == "bv") {
    return std::shared_ptr<Operator>(new Bivar(operator_in));
  } else if (model_type == "bv_matern") {
    return std::shared_ptr<Operator>(new bv_matern(operator_in));
  } else if (model_type == "var1") {
    // VAR(1) bivariate with Cayley reparameterization.
    // Fully implemented in C++; numerical dK via update_all().
    return std::make_shared<RCallback>(operator_in);
  } else {
    throw std::runtime_error("Unknown model.");
  }
};
