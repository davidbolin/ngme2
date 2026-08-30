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

// Probe vectors for the Hutchinson estimators.
//
// Normally these are N_iter Rademacher vectors and the estimators average
// u^T A u over them. When the system is no larger than the probe budget
// (n <= N_iter) that is wasteful and needlessly noisy: n probes along the
// coordinate axes give the trace exactly, for no more work. Scaling them by
// sqrt(n) makes the exact case fall out of the same 1/N_iter averaging the
// stochastic case uses, so no estimator code has to know which mode it is in:
//     sum_i (sqrt(n) e_i)^T A (sqrt(n) e_i) / n = sum_i A_ii = tr(A).
void sparse_llt_solver::ensure_U(unsigned int seed) {
  if (n > 0 && N_iter >= n) {
    if (U_computed && exact_trace && U.rows() == n && U.cols() == n)
      return;
    N_iter = n;
    U = Eigen::MatrixXd::Identity(n, n) * std::sqrt(static_cast<double>(n));
    exact_trace = true;
    U_seed = seed;
    U_computed = true;
    return;
  }
  if (U_computed && !exact_trace && U_seed == seed && U.rows() == n &&
      U.cols() == N_iter)
    return;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  U.resize(n, N_iter);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < N_iter; ++j)
      U(i, j) = dist(rng);
  U = U.unaryExpr(std::ref(myround));
  exact_trace = false;
  U_seed = seed;
  U_computed = true;
}

double sparse_llt_solver::trace(const SparseMatrix<double, 0, int> &M,
                                unsigned int seed) {
  if (QU_computed == 0) {
    ensure_U(seed);
    if (isSymmetric || use_lu) {
      // Symmetric: QU = K^{-1} U. LU mode: solve() is already K^{-1}.
      QU = solve(U);
    } else {
      // For non-symmetric mode, factorization is on Q = K^T K; we need Q^{-1} U
      // (no K^T RHS)
      if (solver_type == 0)
        QU = R_eigen.solve(U);
      else if (solver_type == 1)
        QU = R_ldlt.solve(U);
      else if (solver_type == 2)
        QU = R_supernodal.solve(U);
      else if (solver_type == 3)
        QU = R_cholmod_ldlt.solve(U);
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
  if (isSymmetric || use_lu || K_last.rows() == 0) {
    // QU is already K^{-1} U, so u^T M K^{-1} u estimates tr(K^{-1} M).
    MQU = M * QU;
  } else {
    // Normal equations: QU = (K^T K)^{-1} U, and (K^T K)^{-1} K^T = K^{-1}.
    MQU = (K_last.transpose() * M) * QU;
  }
  double t = 0, t2 = 0;
  for (int i = 0; i < N_iter; i++) {
    double probe = U.col(i).dot(MQU.col(i));
    t += probe;
    t2 += probe * probe;
  }
  double mean = t / N_iter;
  // Spread of the probes: the estimator's own variance is this over N_iter.
  // With the exact (scaled identity) U the trace carries no estimation error,
  // however much the diagonal entries vary, so report zero.
  if (exact_trace || N_iter < 2)
    last_probe_var_ = 0.0;
  else
    last_probe_var_ = (t2 - N_iter * mean * mean) / (N_iter - 1);
  return mean;
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
    ensure_U(seed);
    if (isSymmetric || use_lu) {
      QU = solve(U); // K^{-1} U
    } else {
      // For non-symmetric mode, QU = (K^T K)^{-1} U (no K^T on RHS here)
      if (solver_type == 0)
        QU = R_eigen.solve(U);
      else if (solver_type == 1)
        QU = R_ldlt.solve(U);
      else if (solver_type == 2)
        QU = R_supernodal.solve(U);
      else if (solver_type == 3)
        QU = R_cholmod_ldlt.solve(U);
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
  // In the non-symmetric case, solve() internally applies K^T to the RHS.
  Eigen::MatrixXd A_QU = A * QU;   // n x N_iter
  Eigen::MatrixXd S = solve(A_QU); // Q^{-1} A_eff Q^{-1} U
  Eigen::MatrixXd BS = B * S;      // n x N_iter

  double t = 0.0;
  for (int i = 0; i < N_iter; ++i) {
    t += U.col(i).dot(BS.col(i));
  }
  return t / static_cast<double>(N_iter);
}


// ---------------------------------------------------------------------------
// Selected (Takahashi) inverse.
//
// SimplicialLLT factorizes P Q P^T = L L^T, so Q^{-1} = P^T (L L^T)^{-1} P.
// Writing L = Ltil * diag(L_ii) with Ltil unit-diagonal and D_ii = L_ii^2, the
// entries of B = (L L^T)^{-1} on the pattern of L satisfy, sweeping i upwards
// from the last column,
//     B_ij = delta_ij / D_ii - sum_{k>i, Ltil_ki != 0} Ltil_ki * B_kj.
// Only entries inside pattern(L) are produced; anything outside is genuinely
// unavailable, which is why selinv_trace() reports failure rather than
// silently treating a miss as zero.
// ---------------------------------------------------------------------------
// Factorize privately (once per compute()) so the selected inverse is
// available whatever backend the caller chose.
bool sparse_llt_solver::ensure_selinv_factor() {
  if (!selinv_supported())
    return false;
  // The eigen backend already holds a simplicial factor; reuse it rather than
  // paying for a second factorization.
  if (solver_type == 0)
    return true;
  if (!M_sym_ready)
    return false;
  if (selinv_llt_ready)
    return true;
  selinv_llt.compute(M_sym_);
  if (selinv_llt.info() != Eigen::Success)
    return false;
  selinv_llt_ready = true;
  return true;
}

double sparse_llt_solver::fill_ratio() {
  if (!ensure_selinv_factor() || n <= 0)
    return std::numeric_limits<double>::infinity();
  Eigen::SparseMatrix<double, 0, int> L(
      solver_type == 0 ? Eigen::SparseMatrix<double, 0, int>(R_eigen.matrixL())
                       : Eigen::SparseMatrix<double, 0, int>(selinv_llt.matrixL()));
  return (double)L.nonZeros() / (double)n;
}

bool sparse_llt_solver::build_selinv() {
  if (S_sel_ready)
    return true;
  if (!ensure_selinv_factor())
    return false;

  Eigen::SparseMatrix<double, 0, int> L(
      solver_type == 0 ? Eigen::SparseMatrix<double, 0, int>(R_eigen.matrixL())
                       : Eigen::SparseMatrix<double, 0, int>(selinv_llt.matrixL()));
  Eigen::VectorXd D(n);
  for (int j = 0; j < n; ++j) {
    double d = L.coeff(j, j);
    if (!(d > 0.0))
      return false;
    D(j) = d * d;
  }
  for (int j = 0; j < n; ++j) {
    double d = std::sqrt(D(j));
    for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(L, j); it; ++it)
      it.valueRef() /= d;
  }

  std::vector<std::vector<std::pair<int, double>>> cols(n);
  for (int j = n - 1; j >= 0; --j) {
    std::vector<int> rows;
    for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(L, j); it; ++it)
      if (it.row() > j)
        rows.push_back((int)it.row());
    auto Bget = [&](int a, int b) -> double {
      int hi = a > b ? a : b, lo = a > b ? b : a;
      for (size_t q = 0; q < cols[lo].size(); ++q)
        if (cols[lo][q].first == hi)
          return cols[lo][q].second;
      return 0.0;
    };
    for (int t = (int)rows.size() - 1; t >= 0; --t) {
      int i = rows[t];
      double acc = 0.0;
      for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(L, j); it; ++it) {
        int k = (int)it.row();
        if (k <= j)
          continue;
        acc += it.value() * Bget(i, k);
      }
      cols[j].push_back(std::make_pair(i, -acc));
    }
    double accd = 0.0;
    for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(L, j); it; ++it) {
      int k = (int)it.row();
      if (k <= j)
        continue;
      accd += it.value() * Bget(k, j);
    }
    cols[j].push_back(std::make_pair(j, 1.0 / D(j) - accd));
  }

  std::vector<Eigen::Triplet<double>> trips;
  for (int j = 0; j < n; ++j)
    for (size_t q = 0; q < cols[j].size(); ++q)
      trips.push_back(Eigen::Triplet<double>(cols[j][q].first, j, cols[j][q].second));
  Eigen::SparseMatrix<double, 0, int> S(n, n);
  S.setFromTriplets(trips.begin(), trips.end());
  S.makeCompressed();
  S_sel = S;
  S_sel_ready = true;
  return true;
}

bool sparse_llt_solver::selinv_trace(
    const Eigen::SparseMatrix<double, 0, int> &M, double &out) {
  if (!build_selinv())
    return false;
  // Q^{-1} = P^T B P, so (Q^{-1})_{rc} = B_{P(r),P(c)}.
  const int *pi = (solver_type == 0)
                      ? R_eigen.permutationP().indices().data()
                      : selinv_llt.permutationP().indices().data();
  double acc = 0.0;
  for (int c = 0; c < M.outerSize(); ++c) {
    for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(M, c); it; ++it) {
      int r = (int)it.row();
      int a = pi[r], b = pi[c];
      int hi = std::max(a, b), lo = std::min(a, b);
      double v = S_sel.coeff(hi, lo);
      if (v == 0.0 && hi != lo) {
        // Either a true zero or outside the factor pattern; only the latter is
        // a problem, and distinguishing them costs a pattern probe.
        bool present = false;
        for (Eigen::SparseMatrix<double, 0, int>::InnerIterator si(S_sel, lo); si; ++si)
          if (si.row() == hi) { present = true; break; }
        if (!present)
          return false;
      }
      // tr(Q^{-1} M) = sum_{r,c} (Q^{-1})_{cr} M_{rc}; Q^{-1} is symmetric.
      acc += v * it.value();
    }
  }
  out = acc;
  return true;
}
