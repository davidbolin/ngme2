#include "../include/solver.h"
#include "../include/phase_timing.h"
#include <algorithm>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

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
  { ngme_timing::Scope _s(ngme_timing::rb_qu_solve_us()); QU = solve(U); }
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
  ngme_timing::add(ngme_timing::rb_calls(), 1);
  ngme_timing::Scope _p(ngme_timing::rb_product_us());
  // QU is K^{-1} U, so u^T M K^{-1} u estimates tr(K^{-1} M).
  return reduce_probes(M * QU);
}

Eigen::MatrixXd
sparse_llt_solver::trace_factored_rhs(const SparseMatrix<double, 0, int> &B,
                                      unsigned int seed) {
  ensure_QU(seed);
  return B * QU;
}

double
sparse_llt_solver::trace_factored_with(const SparseMatrix<double, 0, int> &A,
                                       const Eigen::VectorXd &d,
                                       const Eigen::MatrixXd &BQU) {
  ngme_timing::add(ngme_timing::rb_calls(), 1);
  ngme_timing::Scope _p(ngme_timing::rb_product_us());
  Eigen::MatrixXd X = d.asDiagonal() * BQU;
  return reduce_probes(A.transpose() * X);
}

double sparse_llt_solver::trace_factored(const SparseMatrix<double, 0, int> &A,
                                         const Eigen::VectorXd &d,
                                         const SparseMatrix<double, 0, int> &B,
                                         unsigned int seed) {
  ensure_QU(seed);
  ngme_timing::add(ngme_timing::rb_calls(), 1);
  ngme_timing::Scope _p(ngme_timing::rb_product_us());
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
Eigen::MatrixXd
sparse_llt_solver::trace2_lhs(const SparseMatrix<double, 0, int> &A,
                              unsigned int seed) {
  if (QU_computed == 0) {
    ensure_U(seed);
    QU = solve(U); // K^{-1} U
    QU_computed = 1;
  }
  Eigen::MatrixXd A_QU = A * QU;  // A K^{-1} U
  return solve(A_QU);             // K^{-1} A K^{-1} U
}

double sparse_llt_solver::trace2_reduce(const SparseMatrix<double, 0, int> &B,
                                        const Eigen::MatrixXd &S) const {
  Eigen::MatrixXd BS = B * S; // n x N_iter
  double t = 0.0;
  for (int i = 0; i < N_iter; ++i)
    t += U.col(i).dot(BS.col(i));
  return t / static_cast<double>(N_iter);
}

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

  // The factor's own CSC arrays, read in place. Nothing is copied and nothing
  // is rescaled: Ltil_{rj} = L_{rj}/L_{jj} and D_j = L_{jj}^2 are formed on the
  // fly, which removes both a full copy of the factor and a pass over it.
  const Eigen::SparseMatrix<double, 0, int> &L =
      solver_type == 0 ? R_eigen.matrixL().nestedExpression()
                       : selinv_llt->matrixL().nestedExpression();
  if (L.rows() != n || L.cols() != n || !L.isCompressed())
    return false;
  const int *Lp = L.outerIndexPtr();
  const int *Li = L.innerIndexPtr();
  const double *Lx = L.valuePtr();

  // The selected inverse is stored FLAT, sharing L's pattern exactly, so
  // Z(r, c) for r >= c is z[q] at the same offset q that holds L(r, c). That
  // replaces a vector-of-vectors (one heap allocation per column) and lets the
  // final matrix be built by copying L's structure once instead of sorting a
  // triplet array of the whole factor.
  std::vector<double> z(static_cast<size_t>(L.nonZeros()), 0.0);
  // Dense scatter workspaces: pos maps a row of the current column to its slot,
  // lj holds that row's Ltil value. Both turn what used to be a search per
  // access into an array index.
  std::vector<int> pos(n, -1);
  std::vector<double> lj(n, 0.0);
  std::vector<double> acc;

  for (int j = n - 1; j >= 0; --j) {
    const int start = Lp[j], end = Lp[j + 1];
    // Lower-triangular with sorted indices puts the diagonal first. Everything
    // below indexes on that, so check rather than assume.
    if (start >= end || Li[start] != j)
      return false;
    const double Ljj = Lx[start];
    if (!(Ljj > 0.0))
      return false;
    const int m = end - start - 1;

    for (int t = 0; t < m; ++t) {
      const int r = Li[start + 1 + t];
      pos[r] = t;
      lj[r] = Lx[start + 1 + t] / Ljj;
    }
    acc.assign(m, 0.0);

    // acc[a] accumulates -sum_k Ltil_{k,j} Z(r_a, k) over k in this column's
    // rows. Walking column c yields Z(r, c) for r >= c, and that single value
    // serves BOTH the (i = r, k = c) term and the (i = c, k = r) term -- so one
    // pass over the already-computed columns covers every pair, with no search.
    for (int b = 0; b < m; ++b) {
      const int c = Li[start + 1 + b];
      const double lc = lj[c];
      for (int q = Lp[c]; q < Lp[c + 1]; ++q) {
        const int a = pos[Li[q]];
        if (a < 0)
          continue;
        const double v = z[q];
        acc[a] -= lc * v;
        if (Li[q] != c)
          acc[b] -= lj[Li[q]] * v;
      }
    }

    double accd = 0.0;
    for (int t = 0; t < m; ++t) {
      z[start + 1 + t] = acc[t];
      accd += lj[Li[start + 1 + t]] * acc[t];
    }
    z[start] = 1.0 / (Ljj * Ljj) - accd;

    for (int t = 0; t < m; ++t) {
      const int r = Li[start + 1 + t];
      pos[r] = -1;
      lj[r] = 0.0;
    }
  }

  S_sel = L; // pattern; values overwritten below
  std::copy(z.begin(), z.end(), S_sel.valuePtr());
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
