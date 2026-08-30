#include "operator.h"
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

// Unified updater: update K, Z, (optionally) dK, dZ, factorization, and traces
void Operator::update_all(const VectorXd &theta, const UpdateOptions &opts) {
  // 1) Build K and Z once at base theta
  build_KZ(theta);

  // Ensure storage ready
  if ((int)dK.size() != n_theta_K)
    dK.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));
  if ((int)dZ.size() != n_theta_K)
    dZ.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));
  if ((int)d2K.size() != n_theta_K)
    d2K.assign(n_theta_K,
               std::vector<SparseMatrix<double>>(
                   n_theta_K, SparseMatrix<double>(h.size(), h.size())));
  if ((int)d2Z.size() != n_theta_K)
    d2Z.assign(n_theta_K, std::vector<SparseMatrix<double, 0, int>>(
                              n_theta_K, SparseMatrix<double, 0, int>(
                                             h.size(), h.size())));

  // 2) Derivatives: analytic or numeric
  bool have_analytic = false;
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
          dK[j] = (K - K_base) * (1.0 / eps);
        if (opts.compute_dZ)
          dZ[j] = (Z - Z_base) * (1.0 / eps);
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
          dK[j] = (Kp - Km) * (1.0 / (2.0 * eps));
        if (opts.compute_dZ)
          dZ[j] = (Zp - Zm) * (1.0 / (2.0 * eps));
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
  }
  // Storage must be packed before the pattern can be compared (and before the
  // solvers see it).
  if (!K.isCompressed())
    K.makeCompressed();
  // Re-run the symbolic phase exactly when the sparsity pattern moved.
  if (K_pattern_changed() || ngme_counters::cache_disabled()) {
    ngme_counters::bump(ngme_counters::K_analyzes);
    cholK_solver.analyze(K);
    record_K_pattern();
  }
  cholK_solver.compute(K);
  // A failed SparseLU factorization leaves the decomposition unusable: Eigen
  // returns from factorize() before setting up its L store and only guards the
  // solve path with eigen_assert, which is compiled out here, so going on to
  // the traces below would read uninitialised memory and crash the R session.
  // Report it as an ordinary error instead. The Cholesky paths keep their
  // historical behaviour: a failed factorization there is memory-safe, and
  // during SGD it is usually a transient excursion the next iterate recovers
  // from.
  if (cholK_solver.uses_lu() && !cholK_solver.factorization_success()) {
    throw std::runtime_error(
        "Factorization of the operator matrix K failed for model '" +
        generic_type +
        "': K is singular or numerically rank-deficient at the current "
        "parameter values");
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
    for (int j = 0; j < n_theta_K; ++j) {
      for (int k = j; k < n_theta_K; ++k) {
        double t1 = 0.0, t2 = 0.0;
        // -tr(K^{-1} K_k K^{-1} K_j)
        if (dK[k].rows() > 0 && dK[j].rows() > 0) {
          t1 = -cholK_solver.trace2(dK[k], dK[j], opts.trace_seed);
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
