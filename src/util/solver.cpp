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

// QU = K^{-1} U, shared by every estimator below.
//
// solve() is K^{-1} in every mode: the Cholesky of K when K is symmetric, the
// LU of K, or (K^T K)^{-1} K^T when the normal equations are used. So QU is the
// same object however K was factorized, and one formula serves all of them.
// See the note in trace2() for why that matters.
void sparse_llt_solver::ensure_QU(unsigned int seed) {
  if (QU_computed != 0)
    return;
  ensure_U(seed);
  QU = solve(U);
  QU_computed = 1;
}

// Reduce a block MQU = M QU against the probes: mean_i u_i^T (M QU)_i.
double sparse_llt_solver::reduce_probes(const Eigen::MatrixXd &MQU) {
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

double sparse_llt_solver::trace(const SparseMatrix<double, 0, int> &M,
                                unsigned int seed) {
  ensure_QU(seed);
  // QU is K^{-1} U, so u^T M K^{-1} u estimates tr(K^{-1} M).
  return reduce_probes(M * QU);
}

double sparse_llt_solver::trace_factored(const SparseMatrix<double, 0, int> &A,
                                         const Eigen::VectorXd &d,
                                         const SparseMatrix<double, 0, int> &B,
                                         unsigned int seed) {
  ensure_QU(seed);
  // (A^T diag(d) B) QU, right to left, so nothing bigger than n x N_iter is
  // ever built. Mathematically identical to trace(A.transpose() * d.asDiagonal()
  // * B); the association differs, so the two agree to round-off, not bitwise.
  Eigen::MatrixXd X = B * QU;
  X = d.asDiagonal() * X;
  return reduce_probes(A.transpose() * X);
}

// Hutchinson estimator for tr(K^{-1} A K^{-1} B), as
//   QU = K^{-1} U,  S = K^{-1} (A QU),  mean_i u_i^T (B S)_i.
//
// solve() is K^{-1} whichever way K was factorized -- a Cholesky of K when K is
// symmetric, an LU of K, or (K^T K)^{-1} K^T when the normal equations are used
// -- so the same three lines serve every mode and all of them estimate the same
// quantity.
//
// This used to take QU = (K^T K)^{-1} U in the normal-equations mode, which is
// not K^{-1} U: it left an extra K^{-T} in the estimator, so trace2() returned
// tr(K^{-1} A K^{-1} K^{-T} B) rather than the trace asked for, and trace()
// -- which compensated with a K^T A product -- was unbiased but carried around
// a hundred times the standard deviation of this form. Both are why the
// normal-equations mode was previously unusable for a non-symmetric operator.
double sparse_llt_solver::trace2(const SparseMatrix<double, 0, int> &A,
                                 const SparseMatrix<double, 0, int> &B,
                                 unsigned int seed) {
  if (QU_computed == 0) {
    ensure_U(seed);
    QU = solve(U); // K^{-1} U
    QU_computed = 1;
  }

  Eigen::MatrixXd A_QU = A * QU;   // A K^{-1} U
  Eigen::MatrixXd S = solve(A_QU); // K^{-1} A K^{-1} U
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
  if (!selinv_llt)
    selinv_llt.reset(new selinv_llt_t());
  selinv_llt->compute(M_sym_);
  if (selinv_llt->info() != Eigen::Success)
    return false;
  selinv_llt_ready = true;
  return true;
}

double sparse_llt_solver::fill_ratio() {
  if (!ensure_selinv_factor() || n <= 0)
    return std::numeric_limits<double>::infinity();
  // Count through the triangular view's nested expression rather than
  // materialising L: a copy here is the size of the factor itself, and the
  // whole point of this call is to find out whether that size is prohibitive.
  const Eigen::Index nnz =
      solver_type == 0 ? R_eigen.matrixL().nestedExpression().nonZeros()
                       : selinv_llt->matrixL().nestedExpression().nonZeros();
  return (double)nnz / (double)n;
}

bool sparse_llt_solver::build_selinv() {
  if (S_sel_ready)
    return true;
  if (!ensure_selinv_factor())
    return false;

  Eigen::SparseMatrix<double, 0, int> L(
      solver_type == 0 ? Eigen::SparseMatrix<double, 0, int>(R_eigen.matrixL())
                       : Eigen::SparseMatrix<double, 0, int>(selinv_llt->matrixL()));
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
                      : selinv_llt->permutationP().indices().data();
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
