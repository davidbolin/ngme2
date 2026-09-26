#ifndef NGME_ATA_CACHE_H
#define NGME_ATA_CACHE_H

#include "thread_io.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <vector>

// Q = A^T diag(d) A for a sparse A whose SPARSITY PATTERN is fixed while its
// values and d both move.
//
// Eigen redoes the symbolic phase of the triple product on every call. In the
// Gibbs loop A is the operator K, whose pattern is fixed for a whole fit, and d
// is 1/SV, which changes every draw -- so the symbolic work is repeated for
// nothing.
//
// The cache stores, per (row r of A, ordered pair of nonzeros in that row), the
// destination slot in Q's value array. A refill is then one flat scatter:
//
//     for r, for ia in row r, for ib in row r:
//         Q.val[slot[t++]] += d_r * A.val[ia] * A.val[ib]
//
// Deliberately NOT cached: the products A.val[ia] * A.val[ib]. A changes every
// optimizer iteration, so they would need refreshing, and reading them out of A
// is both faster (sequential) and four times leaner than storing them.
//
// A row-major index into A's value array is kept so no transpose of A is needed
// per call; A stays column-major and is read through val_idx_.
class AtDA_cache {
public:
  // A may arrive uncompressed -- the operator's K does. Under uncompressed
  // storage a column's entries run from outer[c] for innerNonZero[c] entries,
  // and the value array has gaps, so nothing may be walked flat.
  static inline int col_begin(const Eigen::SparseMatrix<double, 0, int> &A, int c) {
    return A.outerIndexPtr()[c];
  }
  static inline int col_end(const Eigen::SparseMatrix<double, 0, int> &A, int c) {
    return A.isCompressed() ? A.outerIndexPtr()[c + 1]
                            : A.outerIndexPtr()[c] + A.innerNonZeroPtr()[c];
  }

  // Refill Q with A^T diag(d) A. Returns false if the cache declines the job --
  // pattern too dense for the memory budget -- and the caller must fall back to
  // the direct product. Rebuilds itself whenever A's pattern moves.
  bool refill(const Eigen::SparseMatrix<double, 0, int> &A,
              const Eigen::VectorXd &d,
              Eigen::SparseMatrix<double, 0, int> &Q,
              std::size_t budget_bytes) {
    // Escape hatch: NGME2_NO_ATA falls back to Eigen's triple product without a
    // rebuild. There is no control_opt for this, and the budget check below
    // only guards memory, so this is the only way to take the cache out of a
    // fit if it ever proves wrong. Read once per process, so it cannot be
    // toggled between fits in one session.
    static const bool disabled = std::getenv("NGME2_NO_ATA") != nullptr;
    if (disabled || declined_)
      return false;
    // d indexes rows of A.
    if (d.size() != A.rows())
      return false;
    if (!ready_ || !pattern_matches(A)) {
      if (!build(A, Q, budget_bytes)) {
        declined_ = true;
        return false;
      }
    }
    double *qv = Q.valuePtr();
    std::fill(qv, qv + Q.nonZeros(), 0.0);
    const double *av = A.valuePtr();
    const double *dp = d.data();
    std::size_t t = 0;
    const int n = (int)row_ptr_.size() - 1;
    for (int r = 0; r < n; ++r) {
      const int s0 = row_ptr_[r], s1 = row_ptr_[r + 1];
      if (s0 == s1)
        continue;
      const double dr = dp[r];
      for (int ia = s0; ia < s1; ++ia) {
        const double va = dr * av[val_idx_[ia]];
        for (int ib = s0; ib < s1; ++ib)
          qv[slot_[t++]] += va * av[val_idx_[ib]];
      }
    }
    return true;
  }

private:
  std::vector<int> row_ptr_;   // n+1, row-major offsets into val_idx_
  std::vector<int> val_idx_;   // nnz(A): index into A.valuePtr()
  std::vector<int> col_of_;    // nnz(A): the column each entry belongs to
  std::vector<int> slot_;      // sum_r nnz_r^2: destination in Q.valuePtr()
  std::vector<int> pat_outer_, pat_inner_, pat_nnzcol_;
  int pat_rows_{-1}, pat_cols_{-1};
  bool ready_{false};
  bool declined_{false};

  bool pattern_matches(const Eigen::SparseMatrix<double, 0, int> &A) const {
    if (A.rows() != pat_rows_ || A.cols() != pat_cols_)
      return false;
    if ((int)pat_outer_.size() != A.outerSize() + 1)
      return false;
    if ((int)pat_inner_.size() != A.nonZeros())
      return false;
    if (!std::equal(pat_outer_.begin(), pat_outer_.end(), A.outerIndexPtr()))
      return false;
    if (!A.isCompressed()) {
      if ((int)pat_nnzcol_.size() != A.outerSize())
        return false;
      if (!std::equal(pat_nnzcol_.begin(), pat_nnzcol_.end(), A.innerNonZeroPtr()))
        return false;
    } else if (!pat_nnzcol_.empty()) {
      return false;
    }
    // Row indices in column order, which is the only traversal valid for both
    // storage schemes.
    std::size_t k = 0;
    for (int c = 0; c < A.outerSize(); ++c)
      for (int t = col_begin(A, c); t < col_end(A, c); ++t)
        if (pat_inner_[k++] != A.innerIndexPtr()[t]) return false;
    return true;
  }

  bool build(const Eigen::SparseMatrix<double, 0, int> &A,
             Eigen::SparseMatrix<double, 0, int> &Q,
             std::size_t budget_bytes) {
    const int n = (int)A.rows();
    const int *Ap = A.outerIndexPtr();
    const int *Ai = A.innerIndexPtr();
    const int nnz = (int)A.nonZeros();

    // Row counts, then a counting sort of A's entries into row-major order.
    std::vector<int> cnt(n + 1, 0);
    for (int c = 0; c < A.outerSize(); ++c)
      for (int t = col_begin(A, c); t < col_end(A, c); ++t) cnt[Ai[t] + 1]++;
    // Refuse before allocating if the scatter list would be too large.
    std::size_t pairs = 0;
    for (int r = 0; r < n; ++r) {
      const std::size_t k = (std::size_t)cnt[r + 1];
      pairs += k * k;
    }
    if (pairs * sizeof(int) > budget_bytes)
      return false;

    for (int r = 0; r < n; ++r) cnt[r + 1] += cnt[r];
    row_ptr_.assign(cnt.begin(), cnt.end());
    val_idx_.assign(nnz, 0);
    col_of_.assign(nnz, 0);
    std::vector<int> fill(row_ptr_.begin(), row_ptr_.end() - 1);
    for (int c = 0; c < A.outerSize(); ++c)
      for (int t = col_begin(A, c); t < col_end(A, c); ++t) {
        const int pos = fill[Ai[t]]++;
        val_idx_[pos] = t;
        col_of_[pos] = c;
      }

    // Pattern of Q, built from the contribution list itself so it covers every
    // pair by construction. Deriving it from a product of A instead is wrong:
    // entries that cancel to exactly zero are pruned away and the pair that
    // produced them then has no slot to write to.
    Q.resize(A.cols(), A.cols());
    {
      std::vector<Eigen::Triplet<double>> trips;
      trips.reserve(pairs);
      for (int r = 0; r < n; ++r) {
        const int s0 = row_ptr_[r], s1 = row_ptr_[r + 1];
        for (int ia = s0; ia < s1; ++ia)
          for (int ib = s0; ib < s1; ++ib)
            trips.emplace_back(col_of_[ia], col_of_[ib], 1.0);
      }
      Q.setFromTriplets(trips.begin(), trips.end());
    }
    Q.makeCompressed();

    // Destination slot for every ordered pair, by binary search in Q's column.
    const int *Qp = Q.outerIndexPtr(), *Qi = Q.innerIndexPtr();
    slot_.assign(pairs, 0);
    std::size_t t = 0;
    for (int r = 0; r < n; ++r) {
      const int s0 = row_ptr_[r], s1 = row_ptr_[r + 1];
      for (int ia = s0; ia < s1; ++ia) {
        const int a = col_of_[ia];
        for (int ib = s0; ib < s1; ++ib) {
          const int b = col_of_[ib];
          const int *lo = Qi + Qp[b], *hi = Qi + Qp[b + 1];
          const int *it = std::lower_bound(lo, hi, a);
          if (it == hi || *it != a)
            return false; // pattern of Q does not cover the product
          slot_[t++] = (int)(Qp[b] + (it - lo));
        }
      }
    }
    pat_rows_ = (int)A.rows();
    pat_cols_ = (int)A.cols();
    pat_outer_.assign(A.outerIndexPtr(), A.outerIndexPtr() + A.outerSize() + 1);
    pat_nnzcol_.clear();
    if (!A.isCompressed())
      pat_nnzcol_.assign(A.innerNonZeroPtr(), A.innerNonZeroPtr() + A.outerSize());
    pat_inner_.clear();
    pat_inner_.reserve(A.nonZeros());
    for (int c = 0; c < A.outerSize(); ++c)
      for (int t = col_begin(A, c); t < col_end(A, c); ++t)
        pat_inner_.push_back(A.innerIndexPtr()[t]);
    ready_ = true;
    if (std::getenv("NGME2_ATA_DEBUG")) {
      std::ostringstream msg;
      msg << "[ata] built: A " << A.rows() << 'x' << A.cols()
          << " nnz=" << A.nonZeros() << " pairs=" << slot_.size() << " ("
          << std::fixed << std::setprecision(1)
          << slot_.size() * 4.0 / 1e6 << " MB)\n";
      ngme_io::err() << msg.str();
    }
    return true;
  }
};

#endif
