// Leave-group-out cross-validation: the two C++ entry points R dispatches
// between, plus the kernel the first one is built on.
//
// What both avoid is the shape of cross_validation(), which pays for a
// separate BlockModel construction and a separate Gibbs chain per fold, per
// simulation, per replicate, with every one of those chains refactorizing QQ
// on every sweep.
//
//   group_cv_raw()        the kernel. One full-data Gibbs chain, and each
//                         group handled by a rank-|I| downdate of the factor
//                         sampleW_VY() has already computed -- |I| triangular
//                         solves against an existing factorization rather than
//                         a new one. Returns the leave-group-out conditional
//                         mean and covariance of eta_I per draw, which is
//                         Rao-Blackwellised: R forms the predictive from those
//                         moments and reweights the draws, since they target
//                         p(. | y) rather than p(. | y_-I).
//
//   group_cv_cpp()        runs that kernel over every replicate, OpenMP across
//                         replicates, models built on the main thread because
//                         Rcpp objects must not be touched off it.
//
//   group_cv_exact_cpp()  the other branch entirely: a genuine leave-group-out
//                         chain per fold, with the held-out observations
//                         masked by zeroing their measurement precision. No
//                         downdate and no weights, so nothing to degenerate --
//                         used where the latent field is non-Gaussian and the
//                         weights would. Parallel across folds within a
//                         replicate, each fold restoring the fitted state so
//                         no fold inherits another's and the result does not
//                         depend on scheduling.

// Rinternals.h declares COMPLEX as a function while Accelerate's vDSP.h
// typedefs it, and block.h pulls Accelerate in through the solver. Same guard
// estimate.cpp uses.
#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "block.h"

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

namespace {

// Rows of AZ as (column, value) pairs. AZ depends only on the parameters, which
// are fixed for the whole chain, so this is built once: AZ is column-major, and
// reading a row out of it on every draw would otherwise dominate.
struct RowIndex {
  std::vector<std::vector<std::pair<int, double>>> rows;

  explicit RowIndex(const SparseMatrix<double> &AZ) : rows(AZ.rows()) {
    for (int c = 0; c < AZ.outerSize(); ++c)
      for (SparseMatrix<double>::InnerIterator it(AZ, c); it; ++it)
        rows[it.row()].emplace_back((int)c, it.value());
  }

  void scatter(int r, double *dst) const {
    for (const auto &e : rows[r])
      dst[e.first] = e.second;
  }
  double dot(int r, const double *v) const {
    double acc = 0.0;
    for (const auto &e : rows[r])
      acc += e.second * v[e.first];
    return acc;
  }
};

inline void pack_sym(const MatrixXd &M, double *dst) {
  const int k = (int)M.rows();
  int p = 0;
  for (int j = 0; j < k; ++j)
    for (int i = j; i < k; ++i)
      dst[p++] = M(i, j);
}

} // namespace

// Plain-C++ result, so the chain can run off the main thread: nothing here is
// an R object.
void BlockModel::group_cv_raw(const std::vector<std::vector<int>> &groups,
                              int n, int n_burnin, int chunk_cols,
                              bool allow_inner_threads,
                              std::vector<std::vector<double>> &out_mean,
                              std::vector<std::vector<double>> &out_cov) {
  if (corr_measure)
    throw std::runtime_error(
        "group_cv(): correlated measurement noise is not supported yet. "
        "Removing a group from a correlated likelihood is still a rank-|I| "
        "downdate, but with an effective A_I of A_I + Q_bb^-1 Q_ba A_-I "
        "rather than A_I. Use cross_validation() "
        "for these models.");

  const int n_group = (int)groups.size();
  // A replicate can legitimately hold no groups -- with one fold, only one
  // replicate contains it and the rest get an empty list. Return nothing
  // rather than treating it as an error; the R side already skips them.
  out_mean.clear();
  out_cov.clear();
  if (n_group == 0)
    return;

  std::vector<int> ksize(n_group);
  for (int g = 0; g < n_group; ++g) {
    if (groups[g].empty())
      throw std::runtime_error("group_cv(): group " + std::to_string(g + 1) +
                               " is empty.");
    for (int idx : groups[g])
      if (idx < 0 || idx >= n_obs)
        throw std::runtime_error(
            "group_cv(): observation index out of range in group " +
            std::to_string(g + 1) + ".");
    ksize[g] = (int)groups[g].size();
  }

  out_mean.assign(n_group, {});
  out_cov.assign(n_group, {});
  for (int g = 0; g < n_group; ++g) {
    out_mean[g].assign((size_t)ksize[g] * n, 0.0);
    out_cov[g].assign((size_t)(ksize[g] * (ksize[g] + 1) / 2) * n, 0.0);
  }
  if (n_latent == 0)
    return; // eta_I is identically zero; the zero-filled buffers are the answer

  rao_blackwell = true; // m_f is read off cond_W = QQ^-1 M

  // Every group's observations laid end to end, so one solve serves many
  // groups. Groups are assigned to chunks whole: a group's columns never
  // straddle two solves.
  std::vector<int> col_off(n_group + 1, 0);
  for (int g = 0; g < n_group; ++g)
    col_off[g + 1] = col_off[g] + ksize[g];
  const int total_cols = col_off[n_group];
  if (chunk_cols < 1)
    chunk_cols = 256;
  chunk_cols = std::min(chunk_cols, total_cols);

  std::vector<int> chunk_start;
  {
    int g = 0;
    while (g < n_group) {
      chunk_start.push_back(g);
      int wide = 0;
      while (g < n_group && (wide == 0 || wide + ksize[g] <= chunk_cols)) {
        wide += ksize[g];
        ++g;
      }
    }
    chunk_start.push_back(n_group);
  }
  const int n_chunk = (int)chunk_start.size() - 1;

  if (!all_gaussian)
    burn_in(n_burnin);

  const RowIndex rix(get_AZ());
  const int W_dim = (int)get_AZ().cols();

  for (int it = 0; it < n; ++it) {
    if (!all_gaussian)
      sample_cond_V();
    sampleW_VY();
    // Everything below must run before sample_cond_noise_V(): the downdate has
    // to subtract the measurement precision that actually went into QQ, and
    // that sampler replaces noise_V in place.

    const VectorXd condW = get_cond_W();
    const VectorXd resid = get_residual_part(); // Y - X beta - (noise_V - 1) mu
    const VectorXd s_diag =
        noise_sigma.array().square().matrix().cwiseProduct(noise_V);

    for (int ch = 0; ch < n_chunk; ++ch) {
      const int g0 = chunk_start[ch], g1 = chunk_start[ch + 1];
      const int c0 = col_off[g0], ncol = col_off[g1] - c0;

      // Rebuilt per chunk rather than cached for the whole chain: caching every
      // chunk's dense RHS would cost W_dim * total_cols doubles, which is what
      // chunk_cols exists to bound.
      MatrixXd rhs = MatrixXd::Zero(W_dim, ncol);
      for (int g = g0; g < g1; ++g)
        for (int i = 0; i < ksize[g]; ++i)
          rix.scatter(groups[g][i], rhs.col(col_off[g] - c0 + i).data());

      MatrixXd U = chol_QQ.solve(rhs); // Sigma_f A_I', every column at once

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (allow_inner_threads)
#endif
      for (int g = g0; g < g1; ++g) {
        const int k = ksize[g];
        const int base = col_off[g] - c0;
        MatrixXd Cm(k, k), S(k, k);
        VectorXd mf(k), r(k);
        S.setZero();
        for (int i = 0; i < k; ++i) {
          const int oi = groups[g][i];
          for (int j = 0; j < k; ++j)
            Cm(i, j) = rix.dot(oi, U.col(base + j).data());
          mf(i) = rix.dot(oi, condW.data());
          S(i, i) = s_diag(oi);
          r(i) = resid(oi);
        }
        // S_I - Cm is positive definite by construction: Cm = A_I QQ^-1 A_I'
        // and QQ already carries the +A_I' S_I^-1 A_I term that is being
        // removed, so no fallback is needed here.
        Eigen::LLT<MatrixXd> llt(S - Cm);
        const MatrixXd SCi = llt.solve(MatrixXd::Identity(k, k));
        const VectorXd Sinv_r = S.diagonal().cwiseInverse().cwiseProduct(r);
        const VectorXd m = S * (SCi * (mf - Cm * Sinv_r));
        const MatrixXd Vp = Cm + Cm * SCi * Cm;

        double *mp = out_mean[g].data() + (size_t)k * it;
        for (int i = 0; i < k; ++i)
          mp[i] = m(i);
        pack_sym(Vp, out_cov[g].data() + (size_t)(k * (k + 1) / 2) * it);
      }
    }

    sample_cond_noise_V(true);
  }
}

// Run the whole leave-group-out pass for every replicate in ONE call.
//
// cross_validation() makes n_folds * N_sim * n_replicates separate calls into
// C++, each rebuilding the BlockModel from an R list before it samples
// anything. Here every model is built once, on the main thread (Rcpp objects
// must not be touched off it), the chains then run under OpenMP with
// exceptions captured and re-raised on the main thread, and only the finished
// buffers are copied into R vectors.
//
// [[Rcpp::export]]
Rcpp::List group_cv_cpp(const Rcpp::List &ngme_replicates,
                        const Rcpp::List &groups_per_rep, int n, int n_burnin,
                        unsigned long seed, int num_threads, int chunk_cols) {
  const int n_rep = ngme_replicates.size();
  if (n_rep != groups_per_rep.size())
    Rcpp::stop("group_cv_cpp(): one group list is required per replicate.");
  if (n < 1)
    Rcpp::stop("group_cv_cpp(): n must be at least 1.");

  // 1. Build every model and unpack every group list on the main thread.
  std::vector<std::unique_ptr<BlockModel>> models;
  std::vector<std::vector<std::vector<int>>> groups(n_rep);
  models.reserve(n_rep);
  for (int r = 0; r < n_rep; ++r) {
    models.emplace_back(new BlockModel(Rcpp::as<Rcpp::List>(ngme_replicates[r]),
                                       seed + (unsigned long)r));
    Rcpp::List gl = Rcpp::as<Rcpp::List>(groups_per_rep[r]);
    groups[r].resize(gl.size());
    for (int g = 0; g < gl.size(); ++g) {
      Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(gl[g]);
      groups[r][g].reserve(iv.size());
      for (int i = 0; i < iv.size(); ++i)
        groups[r][g].push_back(iv[i] - 1); // R is 1-based
    }
  }

  // 2. Run the chains. Nothing in here touches the R API.
  std::vector<std::vector<std::vector<double>>> res_mean(n_rep), res_cov(n_rep);
  std::string err;
  const int nt = std::max(1, num_threads);
  // One replicate is the common case, so spend the threads on the groups then;
  // with several replicates the outer axis is the coarser and better one, and
  // nesting the two would oversubscribe.
#ifdef _OPENMP
  const bool par_rep = (n_rep > 1) && (nt > 1);
  const int nt_rep = std::min(nt, n_rep);
#else
  const bool par_rep = false;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt_rep) if (par_rep)
#endif
  for (int r = 0; r < n_rep; ++r) {
    try {
      models[r]->group_cv_raw(groups[r], n, n_burnin, chunk_cols, !par_rep,
                              res_mean[r], res_cov[r]);
    } catch (const std::exception &e) {
#ifdef _OPENMP
#pragma omp critical(ngme_group_cv_exception)
#endif
      if (err.empty())
        err = e.what();
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical(ngme_group_cv_exception)
#endif
      if (err.empty())
        err = "unknown C++ exception";
    }
  }

  if (!err.empty())
    Rcpp::stop("group_cv_cpp(): %s", err.c_str());

  // 3. Copy out, back on the main thread.
  Rcpp::List res(n_rep);
  for (int r = 0; r < n_rep; ++r) {
    const int n_group = (int)res_mean[r].size();
    Rcpp::List m(n_group), v(n_group);
    for (int g = 0; g < n_group; ++g) {
      const int k = (int)(res_mean[r][g].size() / n);
      Rcpp::NumericMatrix mm(k, n);
      std::copy(res_mean[r][g].begin(), res_mean[r][g].end(), mm.begin());
      Rcpp::NumericMatrix vv(k * (k + 1) / 2, n);
      std::copy(res_cov[r][g].begin(), res_cov[r][g].end(), vv.begin());
      m[g] = mm;
      v[g] = vv;
    }
    res[r] = Rcpp::List::create(Rcpp::Named("mean") = m, Rcpp::Named("cov") = v,
                                Rcpp::Named("n_draws") = n);
  }
  return res;
}

// Exact leave-group-out: one Gibbs chain per group, the estimator
// cross_validation() already uses, run entirely inside C++.
//
// [[Rcpp::export]]
Rcpp::List group_cv_exact_cpp(const Rcpp::List &ngme_replicates,
                              const Rcpp::List &groups_per_rep, int n,
                              int n_burnin, unsigned long seed,
                              int num_threads, int n_chains,
                              const Rcpp::List &chain_starts) {
  const int n_rep = ngme_replicates.size();
  if (n_rep != groups_per_rep.size())
    Rcpp::stop("group_cv_exact_cpp(): one group list is required per replicate.");
  if (n < 1)
    Rcpp::stop("group_cv_exact_cpp(): n must be at least 1.");
  if (n_chains < 1)
    Rcpp::stop("group_cv_exact_cpp(): n_chains must be at least 1.");
  const int nt = std::max(1, num_threads);

  // Optional per-chain starting W, one entry per replicate, each a list of
  // vectors (one per chain). Unpacked on the main thread.
  std::vector<std::vector<VectorXd>> starts(n_rep);
  if (chain_starts.size() == n_rep) {
    for (int r = 0; r < n_rep; ++r) {
      Rcpp::List sr = Rcpp::as<Rcpp::List>(chain_starts[r]);
      for (int c = 0; c < sr.size(); ++c)
        starts[r].push_back(Rcpp::as<VectorXd>(sr[c]));
    }
  }

  Rcpp::List res(n_rep);
  for (int r = 0; r < n_rep; ++r) {
    Rcpp::List rep_list = Rcpp::as<Rcpp::List>(ngme_replicates[r]);
    Rcpp::List gl = Rcpp::as<Rcpp::List>(groups_per_rep[r]);
    const int n_group = gl.size();

    std::vector<std::vector<int>> groups(n_group);
    for (int g = 0; g < n_group; ++g) {
      Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(gl[g]);
      groups[g].reserve(iv.size());
      for (int i = 0; i < iv.size(); ++i)
        groups[g].push_back(iv[i] - 1);
    }

    // A replicate holding no folds needs no model at all. max(1, min(nt, 0))
    // is 1, so without this every empty replicate still paid for a full
    // BlockModel construction -- on a fit with many replicates and folds in
    // few of them, that was most of the setup cost.
    if (n_group == 0) {
      res[r] = Rcpp::List::create(Rcpp::Named("eta") = Rcpp::List(0),
                                  Rcpp::Named("n_draws") = n,
                                  Rcpp::Named("n_chains") = n_chains);
      continue;
    }

    // Parallelise over (fold, chain) pairs rather than folds alone. The chains
    // of a fold are independent.
    const int n_task = n_group * n_chains;
    // One model per thread, all built on the main thread.
    const int n_workers = std::max(1, std::min(nt, n_task));
    std::vector<std::unique_ptr<BlockModel>> workers;
    workers.reserve(n_workers);
    for (int t = 0; t < n_workers; ++t) {
      workers.emplace_back(new BlockModel(rep_list, seed + 7919UL * (unsigned long)t));
      // Remember the constructed state -- for a fitted object this is the W and
      // V that estimation left behind, which is far closer to each fold's
      // posterior than the prior. Every fold restores it, so no fold inherits
      // another's state and the result does not depend on scheduling.
      workers.back()->snapshot_state();
    }

    std::vector<std::vector<double>> eta(n_group);
    for (int g = 0; g < n_group; ++g)
      eta[g].assign((size_t)groups[g].size() * n * n_chains, 0.0);
    std::string err;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_workers)               \
    if (n_workers > 1)
#endif
    for (int t = 0; t < n_task; ++t) {
      const int g = t / n_chains, c = t % n_chains;
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      try {
        // Each task seeds from (fold, chain), so the draws do not depend on how
        // the loop was scheduled, and the between-chain spread stays a genuine
        // convergence signal.
        const int k = (int)groups[g].size();
        std::vector<double> one;
        workers[tid]->reseed(seed + 104729UL * (unsigned long)(g + 1) +
                             7907UL * (unsigned long)(c + 1));
        const VectorXd *sw =
            (!starts[r].empty()) ? &starts[r][c % starts[r].size()] : nullptr;
        workers[tid]->loo_chain(groups[g], n, n_burnin, one, sw);
        std::copy(one.begin(), one.end(), eta[g].begin() + (size_t)k * n * c);
      } catch (const std::exception &e) {
#ifdef _OPENMP
#pragma omp critical(ngme_group_cv_exception)
#endif
        if (err.empty()) err = e.what();
      } catch (...) {
#ifdef _OPENMP
#pragma omp critical(ngme_group_cv_exception)
#endif
        if (err.empty()) err = "unknown C++ exception";
      }
    }
    if (!err.empty())
      Rcpp::stop("group_cv_exact_cpp(): %s", err.c_str());

    Rcpp::List out(n_group);
    for (int g = 0; g < n_group; ++g) {
      const int k = (int)groups[g].size();
      // k x (n * n_chains): chains laid end to end, so R can split them.
      Rcpp::NumericMatrix m(k, n * n_chains);
      std::copy(eta[g].begin(), eta[g].end(), m.begin());
      out[g] = m;
    }
    res[r] = Rcpp::List::create(Rcpp::Named("eta") = out,
                                Rcpp::Named("n_draws") = n,
                                Rcpp::Named("n_chains") = n_chains);
  }
  return res;
}
