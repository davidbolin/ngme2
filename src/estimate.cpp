#ifdef USEMKL
#define EIGEN_USE_MKL_ALL
#endif

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include "include/factor_counters.h"
#include "include/timer.h"
#include "ngme.h"
#include "optimizer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <random>
#include <string>

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using std::vector;

using namespace Rcpp;

std::vector<bool> check_conv(const MatrixXd &, const MatrixXd &, int, int,
                             double, double, const std::vector<std::string> &,
                             bool, int, double, const VectorXd &, bool, bool,
                             int n_chains, bool trend_use_tstat,
                             bool use_std_check, double trend_rel_lim,
                             double win_span_iters, bool window_ready,
                             int report_batch,
                             std::vector<double> *rel_drift_out,
                             std::vector<bool> *conv_rhat_out,
                             std::vector<bool> *conv_trend_std_out,
                             std::vector<double> *std_ratio_out,
                             std::vector<double> *slopes_out,
                             std::vector<double> *tstats_out,
                             bool *trend_ready_out);

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List estimate_cpp(const Rcpp::List &R_ngme,
                        const Rcpp::List &control_opt) {
  Eigen::initParallel();
  unsigned long seed = Rcpp::as<unsigned long>(control_opt["seed"]);
  std::mt19937 rng(seed);

  const int iterations = control_opt["iterations"];
  const double max_relative_step = control_opt["max_relative_step"];
  const double max_absolute_step = control_opt["max_absolute_step"];
  const int sampling_strategy = control_opt["sampling_strategy"];
  const bool compute_precond_each_iter = true;
  const int n_repls = R_ngme["n_repls"];
  const bool store_traj = control_opt.containsElementNamed("store_traj")
                              ? Rcpp::as<bool>(control_opt["store_traj"])
                              : true;

  const bool trend_use_tstat =
      control_opt.containsElementNamed("trend_use_tstat")
          ? Rcpp::as<bool>(control_opt["trend_use_tstat"])
          : false;
  const bool use_std_check = control_opt.containsElementNamed("use_std_check")
                                 ? Rcpp::as<bool>(control_opt["use_std_check"])
                                 : true;
  // number of consecutive checkpoints that must pass before stopping
  const int n_conv_batch = control_opt.containsElementNamed("n_conv_batch")
                               ? Rcpp::as<int>(control_opt["n_conv_batch"])
                               : 1;

  Rcpp::List output = R_NilValue;

  vector<vector<VectorXd>> trajs_chains;

  auto timer = std::chrono::steady_clock::now();

  Rcpp::List outputs;
  const bool verbose_enabled = control_opt["verbose"];
  const std::string sgd_method =
      Rcpp::as<std::string>(control_opt["sgd_method"]);

  bool stepsize_decay_enabled = false;
  bool stepsize_decay_on_trend = false;
  int trend_decay_cooldown = 0;
  int stepsize_decay_patience = 0;
  double stepsize_decay_gamma = 1.0;
  double stepsize_decay_min_delta = 0.0;
  int stepsize_decay_warmup = 0;
  if (control_opt.containsElementNamed("stepsize_decay")) {
    std::string decay = Rcpp::as<std::string>(control_opt["stepsize_decay"]);
    stepsize_decay_enabled = (decay == "grad_norm_plateau");
    stepsize_decay_on_trend = (decay == "trend");
  }
  if (stepsize_decay_on_trend) {
    stepsize_decay_gamma =
        Rcpp::as<double>(control_opt["stepsize_decay_gamma"]);
  }
  if (stepsize_decay_enabled) {
    stepsize_decay_patience =
        Rcpp::as<int>(control_opt["stepsize_decay_patience"]);
    stepsize_decay_gamma =
        Rcpp::as<double>(control_opt["stepsize_decay_gamma"]);
    stepsize_decay_min_delta =
        Rcpp::as<double>(control_opt["stepsize_decay_min_delta"]);
    stepsize_decay_warmup = Rcpp::as<int>(control_opt["stepsize_decay_warmup"]);
  }
  if (sgd_method == "bfgs") {
    stepsize_decay_enabled = false;
  }

  const int burnin = control_opt["burnin"];

  int n_min_batch = (control_opt["n_min_batch"]);
  int n_slope_check = (control_opt["n_slope_check"]);
  double std_lim = (control_opt["std_lim"]);
  double trend_lim = (control_opt["trend_lim"]);
  const bool trend_std_conv_check =
      control_opt.containsElementNamed("trend_std_conv_check")
          ? Rcpp::as<bool>(control_opt["trend_std_conv_check"])
          : true;

  int n_chains = (control_opt["n_parallel_chain"]);
  int n_batch = (control_opt["n_batch"]);
  double start_sd = (control_opt["start_sd"]);
  double print_check_info = (control_opt["print_check_info"]);
  double max_R_hat = (control_opt["max_R_hat"]);
  const double trend_rel_lim =
      control_opt.containsElementNamed("trend_rel_lim")
          ? Rcpp::as<double>(control_opt["trend_rel_lim"])
          : 1e-3;
  const bool R_hat_conv_check =
      control_opt.containsElementNamed("R_hat_conv_check")
          ? Rcpp::as<bool>(control_opt["R_hat_conv_check"])
          : true;

  VectorXi num_threads = Rcpp::as<VectorXi>(control_opt["num_threads"]);
  int n_threads_chain = std::max(1, static_cast<int>(num_threads[0]));

#ifdef _OPENMP
  omp_set_max_active_levels(2);
  omp_set_num_threads(num_threads[0] * num_threads[1]);
#else
  (void)n_threads_chain;
#endif

  // When `chain_start` is supplied (one row per chain, in chain order) each
  // chain resumes at its OWN previous endpoint and no jitter is applied.
  MatrixXd chain_start;
  bool has_chain_start = false;
  if (control_opt.containsElementNamed("chain_start") &&
      !Rf_isNull(control_opt["chain_start"])) {
    chain_start = Rcpp::as<MatrixXd>(control_opt["chain_start"]);
    has_chain_start = (chain_start.rows() >= n_chains);
    if (!has_chain_start && chain_start.rows() > 0) {
      Rcpp::warning("chain_start has %d rows but %d chains were requested; "
                    "ignoring it and starting from the fitted parameters.",
                    (int)chain_start.rows(), n_chains);
    }
  }
  // Per-chain LATENT STATE: one model list per chain, differing from R_ngme
  // only in each replicate's W and V. The parameters still come from
  // chain_start; this is what keeps a chain's field and mixing variables its
  // own instead of the chain average.
  Rcpp::List chain_ngme;
  bool has_chain_ngme = false;
  if (control_opt.containsElementNamed("chain_ngme") &&
      !Rf_isNull(control_opt["chain_ngme"])) {
    chain_ngme = Rcpp::as<Rcpp::List>(control_opt["chain_ngme"]);
    has_chain_ngme = (chain_ngme.size() >= n_chains);
  }

  // start_sd exists to disperse chains from a COLD start, so that agreement
  // between them is evidence. Applying it to a continuation destroys exactly
  // the continuity the continuation is for.
  const bool warm_start_no_jitter =
      control_opt.containsElementNamed("warm_start_no_jitter") &&
      Rcpp::as<bool>(control_opt["warm_start_no_jitter"]);

  // init model and optimizer
  vector<std::shared_ptr<Ngme>> ngmes;
  vector<Ngme_optimizer> opt_vec;
  int i = 0;
  for (i = 0; i < n_chains; i++) {
    double sd = (has_chain_start || warm_start_no_jitter || i == 0) ? 0 : start_sd;

    // Not thread-safe using Rcpp::List to init optimizer
    Rcpp::List R_ngme_i =
        has_chain_ngme ? Rcpp::as<Rcpp::List>(chain_ngme[i]) : R_ngme;
    ngmes.push_back(std::make_shared<Ngme>(R_ngme_i, seed + i,
                                           sampling_strategy, num_threads[1],
                                           sd));
    if (has_chain_start) {
      if (chain_start.cols() != ngmes[i]->get_n_params()) {
        Rcpp::stop("chain_start has %d columns but the model has %d "
                   "parameters; the previous fit is not this model.",
                   (int)chain_start.cols(), ngmes[i]->get_n_params());
      }
      ngmes[i]->set_parameter_and_update(chain_start.row(i).transpose(), true);
    }
    opt_vec.push_back(Ngme_optimizer(control_opt, ngmes[i], seed + i));
    if (stepsize_decay_on_trend)
      opt_vec.back().set_stepsize_decay_enabled(true);
    if (verbose_enabled && i > 0) {
      opt_vec.back().set_verbose(false);
    }
  }

  std::vector<std::string> par_names = ngmes[0]->get_par_names();

  // burn in period
  std::atomic<bool> burnin_failed(false);
  std::string burnin_error;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_chain)
#endif
  for (i = 0; i < n_chains; i++) {
    if (burnin_failed.load(std::memory_order_relaxed)) {
      continue;
    }
    try {
      (ngmes[i])->burn_in(burnin);
    } catch (const std::exception &e) {
#ifdef _OPENMP
#pragma omp critical(ngme_parallel_exception)
#endif
      {
        if (!burnin_failed.load(std::memory_order_relaxed)) {
          burnin_failed.store(true, std::memory_order_relaxed);
          burnin_error = e.what();
        }
      }
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical(ngme_parallel_exception)
#endif
      {
        if (!burnin_failed.load(std::memory_order_relaxed)) {
          burnin_failed.store(true, std::memory_order_relaxed);
          burnin_error = "unknown exception";
        }
      }
    }
  }
  if (burnin_failed.load(std::memory_order_relaxed)) {
    Rcpp::stop("C++ exception in parallel burn-in: %s", burnin_error.c_str());
  }

  int n_params = ngmes[0]->get_n_params();
  MatrixXd means(n_batch, n_params);
  MatrixXd vars(n_batch, n_params);

  // for Gelman-Rubin statistic
  MatrixXd batch_sum(n_chains, n_params);
  MatrixXd batch_sq_sum(n_chains, n_params);
  // Non-overlapping sub-batch sums, used to estimate the variance of each
  // chain's batch mean without assuming the iterates are independent.
  int n_sub_batch = 0;
  std::vector<MatrixXd> sub_batch_sum;
  std::vector<int> sub_batch_count;
  VectorXd final_R_hat(n_params);
  final_R_hat.setZero();

  std::vector<Eigen::RowVectorXd> hist_mean;
  std::vector<Eigen::RowVectorXd> hist_var;
  // Per-chain sub-batch means, so R-hat can be computed over the same window
  // the trend test uses instead of over the current batch alone.
  std::vector<MatrixXd> hist_chain_mean;
  // Iteration each history point ends at, so the drift can be expressed as a rate.
  std::vector<int> hist_iter;
  // keep last-diagnostics for reporting
  std::vector<double> last_tstats(n_params, 0.0);
  std::vector<double> last_rel_drift(n_params, 0.0);
  int consec_converged = 0;
  std::vector<bool> last_conv_rhat(n_params, false);
  std::vector<bool> last_conv_trend_std(n_params, false);
  std::vector<double> last_std_ratio(n_params, 0.0);
  std::vector<double> last_slopes(n_params, 0.0);
  bool last_trend_ready = false;
  bool converged_by_param = false;

  std::vector<bool> converge(n_params, false);
  bool all_converge = false;
  int steps = 0;
  int batch_steps = (iterations / n_batch);
  // sqrt(n) blocks of sqrt(n) iterates is the standard batch-means split; it
  // needs enough blocks for a usable variance, so fall back to the naive
  // estimator on very short batches.
  if (batch_steps >= 8) {
    n_sub_batch = std::max(2, (int)std::floor(std::sqrt((double)batch_steps)));
    sub_batch_sum.assign(n_sub_batch, MatrixXd::Zero(n_chains, n_params));
    sub_batch_count.assign(n_sub_batch, 0);
  }

  const int polish_iterations =
      control_opt.containsElementNamed("polish_iterations")
          ? Rcpp::as<int>(control_opt["polish_iterations"])
          : 50;
  const double polish_stepsize_factor =
      control_opt.containsElementNamed("polish_stepsize_factor")
          ? Rcpp::as<double>(control_opt["polish_stepsize_factor"])
          : 0.3;

  int curr_batch = 0;
  double stepsize_decay_scale = 1.0;
  double stepsize_decay_prev_norm = std::numeric_limits<double>::infinity();
  int stepsize_decay_bad_epochs = 0;

  // Enable convergence only when at least one diagnostic is requested.
  const bool param_conv_check = trend_std_conv_check || R_hat_conv_check;
  const bool any_conv_check = param_conv_check;

  while (steps < iterations && !all_converge) {
    MatrixXd mat(n_chains, n_params);
    batch_sum.setZero();
    batch_sq_sum.setZero();
    for (int b = 0; b < n_sub_batch; ++b) {
      sub_batch_sum[b].setZero();
      sub_batch_count[b] = 0;
    }


    // Run batch_steps iterations using unified SGD step (computes grad and,
    // optionally, precond)
    std::atomic<bool> sgd_failed(false);
    std::string sgd_error;
    for (int step = 0; step < batch_steps; step++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_chain)
#endif
      for (i = 0; i < n_chains; i++) {
        if (sgd_failed.load(std::memory_order_relaxed)) {
          continue;
        }
        try {
          // Compute one SGD step; decide whether to compute preconditioner this
          // iter We compute it every iter if compute_precond_each_iter, else
          // let the optimizer refresh at needed cadence
          VectorXd param = opt_vec[i].sgd(0.1, // eps
                                          1,   // one step per loop
                                          max_relative_step, max_absolute_step,
                                          compute_precond_each_iter);
          // if (step == batch_steps - 1) {
#ifdef _OPENMP
#pragma omp critical
#endif
          {
            mat.row(i) = param;
            batch_sum.row(i) += param;
            batch_sq_sum.row(i) += param.array().square().matrix();
            if (n_sub_batch > 0) {
              int sb = (step * n_sub_batch) / batch_steps;
              if (sb >= n_sub_batch)
                sb = n_sub_batch - 1;
              sub_batch_sum[sb].row(i) += param;
              if (i == 0)
                sub_batch_count[sb] += 1;
            }
          }
          // }
        } catch (const std::exception &e) {
#ifdef _OPENMP
#pragma omp critical(ngme_parallel_exception)
#endif
          {
            if (!sgd_failed.load(std::memory_order_relaxed)) {
              sgd_failed.store(true, std::memory_order_relaxed);
              sgd_error = e.what();
            }
          }
        } catch (...) {
#ifdef _OPENMP
#pragma omp critical(ngme_parallel_exception)
#endif
          {
            if (!sgd_failed.load(std::memory_order_relaxed)) {
              sgd_failed.store(true, std::memory_order_relaxed);
              sgd_error = "unknown exception";
            }
          }
        }
      }
      if (sgd_failed.load(std::memory_order_relaxed)) {
        Rcpp::stop("C++ exception in parallel SGD step: %s", sgd_error.c_str());
      }
    }
    steps += batch_steps;

    // compute mean and variance
    means.row(curr_batch) = mat.colwise().mean();
    for (int k = 0; k < n_params; k++) {
      if (n_chains > 1) {
        vars(curr_batch, k) =
            (mat.col(k).array() - means(curr_batch, k)).square().sum() /
            (n_chains - 1);
      } else {
        vars(curr_batch, k) = 0.0;
      }
    }

    if (n_chains > 1) {
      bool batch_converged = false;
      // convergence            // Compute R_hat
      VectorXd R_hat(n_params);
      int n = batch_steps;
      int m = n_chains;

      Eigen::RowVectorXd W(n_params);
      W.setZero();
      Eigen::RowVectorXd B(n_params);
      B.setZero();

      Eigen::RowVectorXd grand_mean(n_params);
      grand_mean.setZero();
      MatrixXd chain_means_mat(n_chains, n_params);

      for (int i = 0; i < n_chains; i++) {
        Eigen::RowVectorXd mean_i = batch_sum.row(i) / n;
        chain_means_mat.row(i) = mean_i;

        Eigen::RowVectorXd s2_i =
            (batch_sq_sum.row(i) - n * mean_i.array().square().matrix()) /
            (n - 1);
        W += s2_i;
        grand_mean += mean_i;
      }
      W /= m;
      grand_mean /= m;

      for (int i = 0; i < n_chains; i++) {
        B += (chain_means_mat.row(i) - grand_mean).array().square().matrix();
      }
      B *= ((double)n / (m - 1));

      // W above is the marginal variance of the iterates. Gelman-Rubin needs
      // the variance the chain mean would have, times n -- i.e. the asymptotic
      // variance sigma^2 * tau. For independent draws the two coincide and
      // R_hat -> 1; for SGD iterates they differ by the integrated
      // autocorrelation time, so using the marginal variance biases R_hat up
      // by sqrt(1 + (tau-1)/n) no matter how well the chains agree.
      //
      // Estimate the asymptotic variance by non-overlapping batch means:
      // split each chain's batch into n_sub_batch blocks, and take
      //     var_hat(chain mean) = var(block means) / n_sub_batch,
      //     W_asym               = n * var_hat(chain mean).
      Eigen::RowVectorXd W_eff = W;
      if (n_sub_batch >= 2) {
        Eigen::RowVectorXd var_of_chain_mean =
            Eigen::RowVectorXd::Zero(n_params);
        int usable_chains = 0;
        for (int i = 0; i < n_chains; i++) {
          Eigen::RowVectorXd blk_mean_sum =
              Eigen::RowVectorXd::Zero(n_params);
          Eigen::RowVectorXd blk_mean_sq =
              Eigen::RowVectorXd::Zero(n_params);
          int n_blk = 0;
          for (int b = 0; b < n_sub_batch; ++b) {
            if (sub_batch_count[b] <= 0)
              continue;
            Eigen::RowVectorXd bm =
                sub_batch_sum[b].row(i) / (double)sub_batch_count[b];
            blk_mean_sum += bm;
            blk_mean_sq += bm.array().square().matrix();
            ++n_blk;
          }
          if (n_blk < 2)
            continue;
          Eigen::RowVectorXd bmean = blk_mean_sum / n_blk;
          Eigen::RowVectorXd s2_blk =
              (blk_mean_sq - n_blk * bmean.array().square().matrix()) /
              (n_blk - 1);
          // variance of the mean of n_blk block means
          var_of_chain_mean += s2_blk / n_blk;
          ++usable_chains;
        }
        if (usable_chains > 0) {
          var_of_chain_mean /= usable_chains;
          W_eff = (double)n * var_of_chain_mean;
        }
      }

      Eigen::RowVectorXd var_hat =
          ((double)(n - 1) / n) * W_eff + (1.0 / n) * B;
      R_hat = (var_hat.array() / W_eff.array()).sqrt().transpose();
      // A degenerate denominator means the chains produced no usable within-
      // chain variation (e.g. a frozen chain). That is not evidence of
      // convergence, so report it as failing rather than as NaN, which every
      // subsequent comparison would silently treat as a pass.
      for (int k = 0; k < n_params; k++) {
        if (!(W_eff(k) > 0.0) || !std::isfinite(R_hat(k)))
          R_hat(k) = std::numeric_limits<double>::infinity();
      }
      final_R_hat = R_hat;

      // Record this batch's sub-batch means (and their across-chain variance)
      // so the trend window can be built from them. With too short a batch to
      // sub-divide, fall back to one point per batch, which is what the test
      // used before.
      if (n_sub_batch >= 2) {
        for (int b = 0; b < n_sub_batch; ++b) {
          if (sub_batch_count[b] <= 0)
            continue;
          Eigen::RowVectorXd m = Eigen::RowVectorXd::Zero(n_params);
          for (int c = 0; c < n_chains; ++c)
            m += sub_batch_sum[b].row(c) / (double)sub_batch_count[b];
          m /= (double)n_chains;
          Eigen::RowVectorXd v = Eigen::RowVectorXd::Zero(n_params);
          for (int c = 0; c < n_chains; ++c) {
            Eigen::RowVectorXd d =
                sub_batch_sum[b].row(c) / (double)sub_batch_count[b] - m;
            v += d.array().square().matrix();
          }
          v /= (double)std::max(1, n_chains - 1);
          MatrixXd cm(n_chains, n_params);
          for (int c = 0; c < n_chains; ++c)
            cm.row(c) = sub_batch_sum[b].row(c) / (double)sub_batch_count[b];
          hist_chain_mean.push_back(cm);
          hist_mean.push_back(m);
          hist_var.push_back(v);
          hist_iter.push_back(steps - batch_steps +
                              (int)std::llround((double)(b + 1) *
                                                (double)batch_steps /
                                                (double)n_sub_batch));
        }
      } else {
        hist_mean.push_back(means.row(curr_batch));
        hist_var.push_back(vars.row(curr_batch));
        hist_chain_mean.push_back(mat);
        hist_iter.push_back(steps);
      }

      // Select the regression window: n_slope_check points, evenly spaced,
      // over the second half of the history.
      const int hist_len = (int)hist_mean.size();
      const bool window_ready = (hist_len >= n_slope_check);
      MatrixXd win_means(n_slope_check, n_params);
      MatrixXd win_vars(n_slope_check, n_params);
      double win_span_iters = 0.0;
      if (window_ready) {
        int lo = hist_len / 2;
        if (hist_len - lo < n_slope_check)
          lo = std::max(0, hist_len - n_slope_check);
        for (int k = 0; k < n_slope_check; ++k) {
          int id = (n_slope_check == 1)
                       ? (hist_len - 1)
                       : lo + (int)std::llround((double)k *
                                                (double)(hist_len - 1 - lo) /
                                                (double)(n_slope_check - 1));
          id = std::min(std::max(id, 0), hist_len - 1);
          win_means.row(k) = hist_mean[id];
          win_vars.row(k) = hist_var[id];
          if (k == 0)
            win_span_iters = -(double)hist_iter[id];
          if (k == n_slope_check - 1)
            win_span_iters += (double)hist_iter[id];
        }
      } else {
        win_means.setZero();
        win_vars.setOnes();
      }


      // Recompute R-hat from the per-chain sub-batch means spanning the second
      // half of the run, rather than from the current batch alone. Same
      // Gelman-Rubin estimator, same batch-means correction for the asymptotic
      // variance -- only the window is longer, which is what makes it stable.
      if (window_ready && (int)hist_chain_mean.size() >= n_slope_check) {
        int L = (int)hist_chain_mean.size();
        int lo = L / 2;
        if (L - lo < n_slope_check)
          lo = std::max(0, L - n_slope_check);
        int n_blk = L - lo;
        int n_win_iter = hist_iter[L - 1] - hist_iter[lo] +
                         (hist_iter.size() > 1
                              ? (hist_iter[1] - hist_iter[0])
                              : 1);
        if (n_blk >= 2 && n_win_iter > 1) {
          MatrixXd chain_mean = MatrixXd::Zero(n_chains, n_params);
          for (int b = lo; b < L; ++b)
            chain_mean += hist_chain_mean[b];
          chain_mean /= (double)n_blk;
          Eigen::RowVectorXd grand = chain_mean.colwise().mean();
          // between-chain: n * variance of the chain means
          Eigen::RowVectorXd Bw = Eigen::RowVectorXd::Zero(n_params);
          for (int c = 0; c < n_chains; ++c)
            Bw += (chain_mean.row(c) - grand).array().square().matrix();
          Bw *= ((double)n_win_iter / (double)(n_chains - 1));
          // within-chain asymptotic variance, by non-overlapping batch means
          Eigen::RowVectorXd W_eff = Eigen::RowVectorXd::Zero(n_params);
          for (int c = 0; c < n_chains; ++c) {
            Eigen::RowVectorXd s2 = Eigen::RowVectorXd::Zero(n_params);
            for (int b = lo; b < L; ++b)
              s2 += (hist_chain_mean[b].row(c) - chain_mean.row(c))
                        .array().square().matrix();
            s2 /= (double)(n_blk - 1);
            W_eff += s2 / (double)n_blk; // var of the chain mean
          }
          W_eff /= (double)n_chains;
          W_eff *= (double)n_win_iter; // -> asymptotic variance
          Eigen::RowVectorXd var_hat =
              ((double)(n_win_iter - 1) / (double)n_win_iter) * W_eff +
              Bw / (double)n_win_iter;
          for (int k = 0; k < n_params; k++) {
            double r = std::sqrt(var_hat(k) / W_eff(k));
            R_hat(k) = (W_eff(k) > 0.0 && std::isfinite(r))
                           ? r
                           : std::numeric_limits<double>::infinity();
          }
          final_R_hat = R_hat;
        }
      }

      if (param_conv_check && curr_batch + 1 >= n_min_batch) {
        converge =
            check_conv(win_means, win_vars, n_slope_check - 1, n_slope_check,
                       std_lim, trend_lim, par_names, print_check_info,
                       batch_steps, max_R_hat, R_hat, trend_std_conv_check,
                       R_hat_conv_check, n_chains, trend_use_tstat,
                       use_std_check, trend_rel_lim, win_span_iters,
                       window_ready, curr_batch + 1, &last_rel_drift,
                       &last_conv_rhat, &last_conv_trend_std, &last_std_ratio,
                       &last_slopes, &last_tstats, &last_trend_ready);
        // COMPACT PROGRESS, for long single runs.
        if (verbose_enabled && !print_check_info) {
          int n_ok = 0;
          int i_rhat = 0, i_t = 0, i_dr = 0;
          for (int i = 0; i < n_params; i++) {
            if (converge[i])
              n_ok++;
            if (final_R_hat(i) > final_R_hat(i_rhat))
              i_rhat = i;
            if (std::abs(last_tstats[i]) > std::abs(last_tstats[i_t]))
              i_t = i;
            if (last_rel_drift[i] > last_rel_drift[i_dr])
              i_dr = i;
          }
          auto nm = [&](int i) {
            return (i < (int)par_names.size()) ? par_names[i] : std::string("?");
          };
          Rcpp::Rcout << "[iter " << steps << "] " << n_ok << "/" << n_params
                      << " converged";
          if (n_ok < n_params) {
            Rcpp::Rcout << "  worst: R_hat " << std::fixed
                        << std::setprecision(3) << final_R_hat(i_rhat) << " ("
                        << nm(i_rhat) << ")";
            if (last_trend_ready) {
              Rcpp::Rcout << ", |t| " << std::setprecision(2)
                          << std::abs(last_tstats[i_t]) << " (" << nm(i_t)
                          << "), drift/100 " << std::setprecision(2)
                          << (100.0 * last_rel_drift[i_dr]) << "% (" << nm(i_dr)
                          << ")";
            }
          }
          Rcpp::Rcout << std::endl;
        }
        bool all_params_ok =
            std::find(begin(converge), end(converge), false) == end(converge);
        // Require the criteria to hold on n_conv_batch consecutive checkpoints.
        // A single passing checkpoint is weak evidence: R-hat bounces around
        // its threshold from batch to batch even long after the chains mix.
        consec_converged = all_params_ok ? consec_converged + 1 : 0;
        batch_converged = (consec_converged >= n_conv_batch);
        if (batch_converged)
          converged_by_param = true;

        // Trend-triggered step-size reduction.  d
        if (stepsize_decay_on_trend && last_trend_ready &&
            trend_decay_cooldown <= 0) {
          bool no_drift = true;
          for (int i = 0; i < n_params; i++)
            if (!(std::abs(last_tstats[i]) <= trend_lim))
              no_drift = false;
          if (no_drift) {
            stepsize_decay_scale *= stepsize_decay_gamma;
            for (int i = 0; i < n_chains; i++)
              opt_vec[i].set_stepsize_decay_scale(stepsize_decay_scale);
            trend_decay_cooldown = n_slope_check; // let a fresh window build
            if (verbose_enabled)
              Rcpp::Rcout << "Trend-triggered stepsize decay at batch "
                          << (curr_batch + 1) << ": scale = "
                          << stepsize_decay_scale << "\n";
          }
        }
        if (trend_decay_cooldown > 0)
          trend_decay_cooldown--;
      }

      // 2. if some parameter converge, stop compute gradient, or slow down the
      // gradient. if (auto_stop)
      //     for (int i=0; i < n_chains; i++) {
      //         ngmes[i]->check_converge(converge);
      //     }
      // Only allow early stop when any diagnostic is active.
      all_converge = any_conv_check && batch_converged;
    }

    if (stepsize_decay_enabled) {
      double mean_grad_norm = 0.0;
      for (int i = 0; i < n_chains; i++) {
        mean_grad_norm += opt_vec[i].get_last_grad_norm();
      }
      mean_grad_norm /= std::max(1, n_chains);

      if (curr_batch < stepsize_decay_warmup) {
        stepsize_decay_prev_norm = mean_grad_norm;
        stepsize_decay_bad_epochs = 0;
      } else {
        if (mean_grad_norm <
            stepsize_decay_prev_norm - stepsize_decay_min_delta) {
          stepsize_decay_bad_epochs = 0;
        } else {
          stepsize_decay_bad_epochs += 1;
        }
        stepsize_decay_prev_norm = mean_grad_norm;

        if (stepsize_decay_bad_epochs >= stepsize_decay_patience) {
          stepsize_decay_bad_epochs = 0;
          stepsize_decay_scale *= stepsize_decay_gamma;
          for (int i = 0; i < n_chains; i++) {
            opt_vec[i].set_stepsize_decay_scale(stepsize_decay_scale);
          }
          if (verbose_enabled) {
            Rcpp::Rcout << "Stepsize decay applied at epoch "
                        << (curr_batch + 1)
                        << ": scale = " << stepsize_decay_scale << "\n";
          }
        }
      }
    }

    curr_batch++;
  }

  // After estimation
  // ****** posterior sampling (sampling each chain..)
  // #pragma omp parallel for schedule(static)
  // for (i=0; i < n_chains; i++) {
  //     ngmes[i]->sampling(100, true);
  // }

  // for testing in the end
  // ngmes[0]->test_in_the_end();

  // ---------------------------------------------------------------------
  // Post-convergence polish.
  //
  // The average is Polyak-Ruppert: the mean of the iterates over the window,
  MatrixXd polish_sum = MatrixXd::Zero(n_chains, n_params);
  int polish_done = 0;
  // steps grows during the polish phase too, so remember where the stopping
  // rule actually fired
  const int steps_at_convergence = steps;
  // Never polish for longer than the fit itself ran. The budget is a fixed
  // count, which is the right shape for a normal fit but disproportionate for a
  // very short one.
  // Also, only polish a run that actually converged.
  const int polish_budget = all_converge ? std::min(polish_iterations, steps) : 0;
  if (polish_budget > 0 && n_chains > 0) {
    double polish_scale = stepsize_decay_scale * polish_stepsize_factor;
    for (i = 0; i < n_chains; i++) {
      opt_vec[i].set_stepsize_decay_enabled(true);
      opt_vec[i].set_stepsize_decay_scale(polish_scale);
    }
    std::atomic<bool> polish_failed(false);
    std::string polish_error;
    for (int step = 0; step < polish_budget; step++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_chain)
#endif
      for (i = 0; i < n_chains; i++) {
        if (polish_failed.load(std::memory_order_relaxed))
          continue;
        try {
          VectorXd param = opt_vec[i].sgd(0.1, 1, max_relative_step,
                                          max_absolute_step,
                                          compute_precond_each_iter);
#ifdef _OPENMP
#pragma omp critical
#endif
          { polish_sum.row(i) += param; }
        } catch (const std::exception &e) {
#ifdef _OPENMP
#pragma omp critical(ngme_parallel_exception)
#endif
          {
            if (!polish_failed.load(std::memory_order_relaxed)) {
              polish_failed.store(true, std::memory_order_relaxed);
              polish_error = e.what();
            }
          }
        }
      }
      if (polish_failed.load(std::memory_order_relaxed))
        break;
      polish_done++;
      steps++;
    }
    if (polish_failed.load(std::memory_order_relaxed))
      Rcpp::warning("ngme: post-convergence polish stopped early: %s",
                    polish_error.c_str());
    if (verbose_enabled && polish_done > 0)
      Rcpp::Rcout << "Polish phase: " << polish_done
                  << " iterations at stepsize scale " << polish_scale << "\n";
  }

  // generate outputs
  // final_params keeps each chain's own parameter vector, which is what a later
  // call needs to resume that chain
  MatrixXd final_params(n_chains, n_params);
  for (i = 0; i < n_chains; i++) {
    // Average of the iterates: over the polish window when there was one,
    // otherwise over the last batch as before.
    Eigen::VectorXd avg_param =
        (polish_done > 0)
            ? (polish_sum.row(i) / polish_done).transpose()
            : (batch_sum.row(i) / batch_steps).transpose();
    ngmes[i]->set_parameter_and_update(avg_param, false);
    final_params.row(i) = avg_param.transpose();

    outputs.push_back(ngmes[i]->output());
    if (store_traj) {
      opt_vec[i].record_current_state();
      trajs_chains.push_back(opt_vec[i].get_trajs());
    }
  }
  if (all_converge) {
    Rcpp::Rcout << "Reach convergence in " << steps_at_convergence
                << " iterations." << std::endl;
    if (polish_done > 0)
      Rcpp::Rcout << "Polish: " << polish_done
                  << " further iterations averaged at stepsize scale "
                  << (stepsize_decay_scale * polish_stepsize_factor) << ".\n";
  }

  // A run that exhausts its budget without converging used to end silently,
  // which is the dangerous failure mode: the estimates look like any other
  // result. Say so, and name the parameters that were still failing.
  const bool warn_no_convergence =
      control_opt.containsElementNamed("warn_no_convergence")
          ? Rcpp::as<bool>(control_opt["warn_no_convergence"])
          : true;
  if (!all_converge && any_conv_check && n_chains > 1 && warn_no_convergence) {
    std::string failing;
    for (int i = 0; i < n_params && i < (int)par_names.size(); i++) {
      bool ok = true;
      if (R_hat_conv_check && !last_conv_rhat[i])
        ok = false;
      if (trend_std_conv_check && (!last_trend_ready || !last_conv_trend_std[i]))
        ok = false;
      if (!ok) {
        if (!failing.empty())
          failing += ", ";
        failing += par_names[i];
      }
    }
    if (failing.empty())
      failing = "(criteria never held on " + std::to_string(n_conv_batch) +
                " consecutive checkpoints)";
    Rcpp::warning("ngme: convergence was NOT reached in %d iterations. Still "
                  "failing: %s. Estimates may be unreliable -- increase "
                  "'iterations', or inspect traceplot().",
                  steps_at_convergence, failing.c_str());
  }

  if (all_converge && n_chains > 1) {
    Rcpp::Rcout << "Convergence criteria summary:\n";

    {
      if (R_hat_conv_check) {
        Rcpp::Rcout << "  - R_hat threshold: max_R_hat = " << max_R_hat << "\n";
      }
      if (trend_std_conv_check) {
        Rcpp::Rcout << "  - Trend/Std thresholds: std_lim = " << std_lim
                    << ", trend_lim = " << trend_lim
                    << ", window (stop points) = " << n_slope_check << "\n";
      }

      Rcpp::Rcout << "  Per-parameter status (R_hat | std/mean | slope):\n";
      for (int i = 0; i < n_params; i++) {
        Rcpp::Rcout << "    * " << par_names[i] << ": ";
        bool printed_any = false;
        if (R_hat_conv_check) {
          Rcpp::Rcout << "R_hat=" << std::fixed << std::setprecision(3)
                      << final_R_hat(i)
                      << (last_conv_rhat[i] ? " (ok)" : " (fail)");
          printed_any = true;
        }
        if (trend_std_conv_check && last_trend_ready) {
          if (printed_any)
            Rcpp::Rcout << "; ";
          Rcpp::Rcout << "std/mean=" << std::setprecision(3)
                      << last_std_ratio[i]
                      << (last_conv_trend_std[i] ? " (ok)" : " (fail)") << ", ";
          Rcpp::Rcout << "t=" << std::setprecision(3) << last_tstats[i] << ", ";
        Rcpp::Rcout << "slope=" << std::setprecision(3) << last_slopes[i]
                    << (std::abs(last_slopes[i]) <= trend_lim ? " (ok)"
                                                              : " (fail)");
        }
        Rcpp::Rcout << "\n";
      }
    }
  }

  if (store_traj) {
    outputs.attr("opt_traj") = trajs_chains;
  } else {
    outputs.attr("opt_traj") = R_NilValue;
  }
  if (n_chains > 1)
    outputs.attr("R_hat") = final_R_hat;
  outputs.attr("chain_params") = Rcpp::wrap(final_params);
  outputs.attr("par_names") = Rcpp::wrap(par_names);
  return outputs;
}

// A is the horizontal concat version of latent A
// Rcpp::List sampling_cpp_given_A(const Rcpp::List& ngme_replicate, int n,
// bool posterior, unsigned long seed, const Eigen::SparseMatrix<double>& A) {
//     std::mt19937 rng (seed);
//     BlockModel block (ngme_replicate, rng());
//     return block.sampling(n, posterior, A);
// }

// [[Rcpp::export]]
Rcpp::List sampling_cpp(const Rcpp::List &ngme_replicate, int n, int n_burnin,
                        bool posterior, unsigned long seed) {
  std::mt19937 rng(seed);
  BlockModel block(ngme_replicate, seed);

  return block.sampling(n, n_burnin, posterior);
}

/*
    For checking convergence of parallel chains
    data is n_iters () * n_params (how many params in total)
*/
std::vector<bool>
check_conv(const MatrixXd &means, const MatrixXd &vars, int curr_batch,
           int n_slope_check, double std_lim, double trend_lim,
           const std::vector<std::string> &par_names, bool print_check_info,
           int batch_steps, double max_R_hat, const VectorXd &R_hat,
           bool trend_std_conv_check, bool R_hat_conv_check, int n_chains,
           bool trend_use_tstat, bool use_std_check, double trend_rel_lim,
           double win_span_iters, bool window_ready, int report_batch,
           std::vector<double> *rel_drift_out, std::vector<bool> *conv_rhat_out,
           std::vector<bool> *conv_trend_std_out,
           std::vector<double> *std_ratio_out, std::vector<double> *slopes_out,
           std::vector<double> *tstats_out, bool *trend_ready_out) {
  int n_params = means.cols();
  std::vector<bool> conv(n_params, false);

  std::vector<bool> conv_rhat(n_params, true);
  std::vector<bool> conv_trend_std(n_params, true);
  std::vector<double> std_ratio(n_params, 0.0);
  std::vector<double> slopes(n_params, 0.0);
  std::vector<double> tstats(n_params, 0.0);
  std::vector<double> rel_drift(n_params, 0.0);

  bool trend_ready = window_ready && (curr_batch + 1 >= n_slope_check);

  (void)batch_steps; // currently unused but kept for signature stability

  if (!trend_std_conv_check && !R_hat_conv_check) {
    std::fill(conv.begin(), conv.end(), true);
    if (conv_rhat_out)
      *conv_rhat_out = conv_rhat;
    if (conv_trend_std_out)
      *conv_trend_std_out = conv_trend_std;
    if (std_ratio_out)
      *std_ratio_out = std_ratio;
    if (slopes_out)
      *slopes_out = slopes;
    if (tstats_out)
      *tstats_out = tstats;
    if (trend_ready_out)
      *trend_ready_out = trend_ready;
    return conv;
  }

  if (R_hat_conv_check) {
    for (int i = 0; i < n_params; i++) {
      conv_rhat[i] = (R_hat(i) <= max_R_hat);
    }
  }

  if (trend_std_conv_check && trend_ready) {
    for (int i = 0; i < n_params; i++) {
      std_ratio[i] = std::sqrt(vars(curr_batch, i)) /
                     (std::abs(means(curr_batch, i)) + 1e-5);
      // The raw coefficient of variation is not a convergence statistic: it
      // conflates "the chains disagree" with "this parameter is imprecisely
      // determined", so a weakly identified parameter (nu) can never pass no
      // matter how well mixed the chains are. Off by default; R-hat already
      // measures between- vs within-chain disagreement in a scale-free way.
      if (use_std_check && std_ratio[i] > std_lim) {
        conv_trend_std[i] = false;
      }
    }

    // Build design matrix once
    MatrixXd B(n_slope_check, 2);
    B.col(0) = VectorXd::Ones(n_slope_check);
    for (int i = 0; i < n_slope_check; i++) {
      B(i, 1) = i;
    }

    for (int i = 0; i < n_params; i++) {
      VectorXd mean =
          means.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1);
      VectorXd Sigma_inv =
          vars.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1)
              .cwiseInverse();

      MatrixXd Q = B.transpose() * Sigma_inv.asDiagonal() * B;
      Vector2d beta =
          Q.llt().solve(B.transpose() * Sigma_inv.asDiagonal() * mean);

      slopes[i] = beta(1);

      // means(b, i) is the average over chains, so its sampling variance is
      // vars(b, i) / n_chains; Sigma_inv above omits that factor. It cancels in
      // beta but not in its covariance, so put it back here.
      // cov(beta) = Q^{-1} with the corrected weights => scale Q by n_chains.
      Matrix2d Qs = Q * (double)n_chains;
      double var_slope = Qs.inverse()(1, 1);
      tstats[i] = (var_slope > 0 && std::isfinite(var_slope))
                      ? beta(1) / std::sqrt(var_slope)
                      : std::numeric_limits<double>::infinity();

      // Scale-free drift test: is the slope distinguishable from zero relative
      // to its own standard error?
      double stat = trend_use_tstat ? std::abs(tstats[i]) : std::abs(beta(1));

      // Also require the drift to be material as a rate per 100 iterations.
      double span_pts = (double)(n_slope_check - 1);
      double scale = std::abs(means(curr_batch, i)) +
                     std::sqrt(std::max(0.0, vars(curr_batch, i))) + 1e-10;
      double per_iter = (win_span_iters > 0.0)
                            ? (std::abs(beta(1)) * span_pts / win_span_iters)
                            : 0.0;
      rel_drift[i] = 100.0 * per_iter / scale; // fraction of scale per 100 iters
      if (stat > trend_lim && rel_drift[i] > trend_rel_lim) {
        conv_trend_std[i] = false;
      }
    }
  }

  // Every enabled diagnostic must pass
  for (int i = 0; i < n_params; i++) {
    bool passed = true;
    if (trend_std_conv_check)
      passed = passed && trend_ready && conv_trend_std[i];
    if (R_hat_conv_check)
      passed = passed && conv_rhat[i];
    conv[i] = passed;
  }

  // copy out results for reporting
  if (conv_rhat_out)
    *conv_rhat_out = conv_rhat;
  if (conv_trend_std_out)
    *conv_trend_std_out = conv_trend_std;
  if (std_ratio_out)
    *std_ratio_out = std_ratio;
  if (slopes_out)
    *slopes_out = slopes;
  if (tstats_out)
    *tstats_out = tstats;
  if (rel_drift_out)
    *rel_drift_out = rel_drift;
  if (trend_ready_out)
    *trend_ready_out = trend_ready;

  if (print_check_info) {
    // curr_batch indexes the regression WINDOW (the caller passes the selected
    // points, not the whole history), so the checkpoint number for the reader
    // has to be handed in separately. Otherwise every checkpoint prints as
    // "stop n_slope_check".
    Rcpp::Rcout << "\nstop " << report_batch << ":\n";

    const int label_width = 11; // width for the row label (e.g., "R_hat:")
    const int col_width = 10;   // width for each value/parameter name
    int line_width = label_width + n_params * (col_width + 1);

    auto print_separator = [&]() {
      Rcpp::Rcout << std::string(line_width, '-') << "\n";
    };

    print_separator();

    Rcpp::Rcout << std::setw(label_width) << std::left << "Param:";
    for (const auto &name : par_names) {
      Rcpp::Rcout << " " << std::setw(col_width) << std::left << name;
    }
    Rcpp::Rcout << "\n";

    print_separator();

    if (R_hat_conv_check) {
      Rcpp::Rcout << std::setw(label_width) << std::left << "R_hat:";
      for (int i = 0; i < n_params; i++) {
        Rcpp::Rcout << " " << std::setw(col_width) << std::fixed
                    << std::setprecision(3) << std::left << R_hat(i);
      }
      Rcpp::Rcout << "\n";
    }

    if (trend_std_conv_check && trend_ready) {
      Rcpp::Rcout << std::setw(label_width) << std::left << "std/mean:";
      for (int i = 0; i < n_params; i++) {
        Rcpp::Rcout << " " << std::setw(col_width) << std::fixed
                    << std::setprecision(3) << std::left << std_ratio[i];
      }
      Rcpp::Rcout << "\n";

      Rcpp::Rcout << std::setw(label_width) << std::left << "trend:";
      for (int i = 0; i < n_params; i++) {
        Rcpp::Rcout << " " << std::setw(col_width) << std::fixed
                    << std::setprecision(3) << std::left << slopes[i];
      }
      Rcpp::Rcout << "\n";

      Rcpp::Rcout << std::setw(label_width) << std::left << "drift/100:";
      for (int i = 0; i < n_params; i++) {
        Rcpp::Rcout << " " << std::setw(col_width) << std::scientific
                    << std::setprecision(2) << std::left << rel_drift[i];
      }
      Rcpp::Rcout << std::fixed << "\n";

      Rcpp::Rcout << std::setw(label_width) << std::left << "t_slope:";
      for (int i = 0; i < n_params; i++) {
        Rcpp::Rcout << " " << std::setw(col_width) << std::fixed
                    << std::setprecision(3) << std::left << tstats[i];
      }
      Rcpp::Rcout << "\n";
    }

    print_separator();
    Rcpp::Rcout << "\n";
  }

  return conv;
}

// check convergence of parallel chains
// std.lim, trend.lim
//   if(!is.null(dim(m))){
//     n.test <- dim(m)[2]
//     N <- dim(m)[1]
//     output <- rep(FALSE,n.test)
//     if(N>3){
//       n.points <- min(N,4)
//       B <- cbind(rep(1,n.points),0:(n.points-1))
//       for(i in 1:n.test){
//         std.satisfied <- sqrt(sigma2[N,i])/abs(m[N,1])<std.lim
//         Sigma <- diag(sigma2[(N-n.points+1):N,i])
//         Q <- solve(t(B)%*%solve(Sigma,B))
//         beta <- Q%*%(t(B)%*%solve(Sigma,m[(N-n.points+1):N,i]))
//         slope.satisfied <-
//         abs(beta[2])-2*sqrt(Q[2,2])<trend.lim*abs(beta[1]) #no significant
//         trend output[i] = std.satisfied&slope.satisfied
//       }
//     }
//     return(output)
// [[Rcpp::export]]
int get_openmp_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 0;
#endif
}

// [[Rcpp::export]]
bool has_pardiso() {
#ifdef USEMKL
  return true;
#else
  return false;
#endif
}

// Read (and optionally reset) the factorization counters. Used by the test
// suite to check that QQ is assembled/factorized once per optimizer iteration
// instead of once per Gibbs draw, and that the symbolic phase only reruns when
// a sparsity pattern actually changes.
// [[Rcpp::export]]
Rcpp::List ngme_factor_counters(bool reset = false) {
  Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("QQ_builds") = static_cast<double>(
          ngme_counters::QQ_builds.load(std::memory_order_relaxed)),
      Rcpp::Named("QQ_analyzes") = static_cast<double>(
          ngme_counters::QQ_analyzes.load(std::memory_order_relaxed)),
      Rcpp::Named("K_analyzes") = static_cast<double>(
          ngme_counters::K_analyzes.load(std::memory_order_relaxed)));
  if (reset)
    ngme_counters::reset_all();
  return out;
}
