#ifdef USEMKL
#define EIGEN_USE_MKL_ALL
#endif

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include "include/factor_counters.h"
#include "include/thread_io.h"
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
#include <deque>
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
                             int n_chains,
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
  const bool trend_std_conv_check_req =
      control_opt.containsElementNamed("trend_std_conv_check")
          ? Rcpp::as<bool>(control_opt["trend_std_conv_check"])
          : true;

  int n_chains = (control_opt["n_parallel_chain"]);
  int n_batch = (control_opt["n_batch"]);
  double start_sd = (control_opt["start_sd"]);
  double print_check_info = (control_opt["print_check_info"]);
  double max_R_hat = (control_opt["max_R_hat"]);
  // "stationarity" stops when the chains have settled, leaving the remaining
  // Monte Carlo error to the Polyak-Ruppert average of the iterates rather than
  // trying to grind it down by running longer. "drift" is the older rule, which
  // required the fitted drift rate itself to fall below trend_rel_lim.
  const std::string conv_criterion =
      control_opt.containsElementNamed("conv_criterion")
          ? Rcpp::as<std::string>(control_opt["conv_criterion"])
          : std::string("stationarity");
  const bool stationarity_mode = (conv_criterion == "stationarity");
  // In stationarity mode the drift-rate clause is replaced, not added to: it is
  // the clause whose cost this mode exists to avoid.
  const bool trend_std_conv_check =
      trend_std_conv_check_req && !stationarity_mode;
  // Length of the stationarity comparison window, in sub-batches.
  // 0 = derive it from the measured correlation length of the chains.
  const int stationarity_window =
      control_opt.containsElementNamed("stationarity_window")
          ? Rcpp::as<int>(control_opt["stationarity_window"])
          : 0;
  const double stationarity_eff =
      control_opt.containsElementNamed("stationarity_eff")
          ? Rcpp::as<double>(control_opt["stationarity_eff"])
          : 8.0;
  // MAGNITUDE, not significance (2026-09-03). See the stationarity test below.
  const double polish_floor_frac =
      control_opt.containsElementNamed("polish_floor_frac")
          ? Rcpp::as<double>(control_opt["polish_floor_frac"])
          : 0.25;
  const double stationarity_ratio_lim =
      control_opt.containsElementNamed("stationarity_ratio_lim")
          ? Rcpp::as<double>(control_opt["stationarity_ratio_lim"])
          : 1.0;
  const double stationarity_dir_lim =
      control_opt.containsElementNamed("stationarity_dir_lim")
          ? Rcpp::as<double>(control_opt["stationarity_dir_lim"])
          : 0.5;
  const int stationarity_min_checks =
      control_opt.containsElementNamed("stationarity_min_checks")
          ? Rcpp::as<int>(control_opt["stationarity_min_checks"])
          : 6;
  const double schedule_arm_z =
      control_opt.containsElementNamed("schedule_arm_z")
          ? Rcpp::as<double>(control_opt["schedule_arm_z"])
          : 2.0;
  const bool schedule_auto_start =
      control_opt.containsElementNamed("schedule_auto_start")
          ? Rcpp::as<bool>(control_opt["schedule_auto_start"])
          : false;
  const bool mc_se_conv_check =
      control_opt.containsElementNamed("mc_se_conv_check")
          ? Rcpp::as<bool>(control_opt["mc_se_conv_check"])
          : true;
  const double mc_se_lim = control_opt.containsElementNamed("mc_se_lim")
                               ? Rcpp::as<double>(control_opt["mc_se_lim"])
                               : 0.5;
  const int n_settle_checks =
      control_opt.containsElementNamed("n_settle_checks")
          ? Rcpp::as<int>(control_opt["n_settle_checks"])
          : 3;
  const int max_stepsize_decays =
      control_opt.containsElementNamed("max_stepsize_decays")
          ? Rcpp::as<int>(control_opt["max_stepsize_decays"])
          : 6;
  const double stepsize_decay_precision_gamma =
      control_opt.containsElementNamed("stepsize_decay_precision_gamma")
          ? Rcpp::as<double>(control_opt["stepsize_decay_precision_gamma"])
          : 0.5;
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
  // chain resumes at its own previous endpoint and no jitter is applied.
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
  // Per-chain latent state: one model list per chain, differing from R_ngme
  // only in each replicate's W and V. The parameters still come from
  // chain_start.
  Rcpp::List chain_ngme;
  bool has_chain_ngme = false;
  if (control_opt.containsElementNamed("chain_ngme") &&
      !Rf_isNull(control_opt["chain_ngme"])) {
    chain_ngme = Rcpp::as<Rcpp::List>(control_opt["chain_ngme"]);
    has_chain_ngme = (chain_ngme.size() >= n_chains);
  }

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
    // With auto-arming the schedule must be inert until the transient is over,
    // not merely re-based when it fires. poly_decay() defaults burnin_iter to 0,
    // so without this the decay runs from iteration 0, freezes the chains short
    // of the optimum, and the drift can then never fall far enough to arm.
    if (schedule_auto_start)
      opt_vec.back().set_schedule_burnin(std::numeric_limits<int>::max() / 2);
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
  std::vector<std::vector<double>> conv_diag; // see the checkpoint block below
  VectorXd se_stat = VectorXd::Constant(n_params, NA_REAL);
  VectorXd se_smooth = VectorXd::Constant(n_params, NA_REAL);
  VectorXd mc_se = VectorXd::Constant(n_params, NA_REAL);
  // Smoothed mc_se.  A single checkpoint's batch-means estimate swings by an
  // order of magnitude between adjacent checkpoints so gating on the raw value
  // lets two consecutive lucky readings declare convergence.
  VectorXd mc_smooth = VectorXd::Constant(n_params, NA_REAL);
  int n_stepsize_decays = 0;
  // Index into the sub-batch history below which points are stale.  Changing the
  // step size changes the stationary distribution the chains live in, so every
  // iterate taken before the change belongs to a different distribution and must
  // not enter any window.
  int hist_floor = 0;
  // Direction history for the flat-direction test. Must live outside the
  // per-checkpoint block.
  std::vector<std::deque<double>> stat_dir_hist;
  std::vector<bool> stat_flat;
  bool stat_dir_init = false;
  // The step-size schedule is armed when the transient ends, not at a preset
  // iteration. The trend test asks exactly whether the parameters have stopped
  // moving systematically. A run that is settled but noisy arms the schedule
  // and gets its noise damped, where a rule gated on R-hat would deadlock.
  bool schedule_armed = false;
  int schedule_ready_count = 0;
  // Checkpoints still to wait through after a step-size change.  The window is
  // rebuilt from scratch then, so the first checkpoints after a decay measure
  // mc_se from the minimum number of blocks and are at their noisiest.
  // Waiting lets the new step size accumulate enough blocks to be judged on.
  int settle_left = 0;
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
  // sqrt(n) blocks of sqrt(n) iterates is the standard batch-means split. It
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
  // The polish decays its step size rather than dropping it once.
  const double polish_decay_alpha =
      control_opt.containsElementNamed("polish_decay_alpha")
          ? Rcpp::as<double>(control_opt["polish_decay_alpha"])
          : 0.6;
  const double polish_decay_t0 =
      control_opt.containsElementNamed("polish_decay_t0")
          ? Rcpp::as<double>(control_opt["polish_decay_t0"])
          : 100.0;
  const double polish_stepsize_factor =
      control_opt.containsElementNamed("polish_stepsize_factor")
          ? Rcpp::as<double>(control_opt["polish_stepsize_factor"])
          : 0.3;

  int curr_batch = 0;
  double stepsize_decay_scale = 1.0;
  double stepsize_decay_prev_norm = std::numeric_limits<double>::infinity();
  int stepsize_decay_bad_epochs = 0;

  // Enable convergence only when at least one diagnostic is requested.
  const bool param_conv_check =
      trend_std_conv_check || R_hat_conv_check || stationarity_mode;
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
          // Compute one SGD step and decide whether to compute preconditioner this
          // iter. We compute it every iter if compute_precond_each_iter, else
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

      // Gelman-Rubin needs the asymptotic variance, not the marginal variance
      // of the iterates: for autocorrelated draws the two differ by tau and
      // R_hat is biased up. Estimate it by non-overlapping batch means,
      //     W_asym = n * var(block means) / n_sub_batch.
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
      // sub-divide, fall back to one point per batch.
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
      // enough FRESH points (post any step-size change) to fill the window
      const bool window_ready = (hist_len - hist_floor >= n_slope_check);
      MatrixXd win_means(n_slope_check, n_params);
      MatrixXd win_vars(n_slope_check, n_params);
      double win_span_iters = 0.0;
      if (window_ready) {
        int lo = std::max(hist_floor, (hist_len + hist_floor) / 2);
        if (hist_len - lo < n_slope_check)
          lo = std::max(hist_floor, hist_len - n_slope_check);
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
      // half of the run, rather than from the current batch alone.
      if (window_ready && (int)hist_chain_mean.size() - hist_floor >= n_slope_check) {
        int L = (int)hist_chain_mean.size();
        int lo = std::max(hist_floor, (L + hist_floor) / 2);
        if (L - lo < n_slope_check)
          lo = std::max(hist_floor, L - n_slope_check);
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


      {
        se_stat = VectorXd::Constant(n_params, NA_REAL);
        if (sgd_method == "precond_sgd") {
          MatrixXd P = MatrixXd::Zero(n_params, n_params);
          for (int c = 0; c < n_chains; c++)
            P += opt_vec[c].get_preconditioner();
          P /= (double)n_chains;
          Eigen::FullPivLU<MatrixXd> lu(P);
          if (lu.isInvertible()) {
            MatrixXd Pinv = lu.inverse();
            for (int k = 0; k < n_params; k++)
              se_stat(k) = (Pinv(k, k) > 0) ? std::sqrt(Pinv(k, k)) : NA_REAL;
          }
        }
        // se_stat is from the complete-data information, so it understates the
        // true standard error by the missing-information fraction.
        int Lh = (int)hist_chain_mean.size();
        int loh = std::max(hist_floor, (Lh + hist_floor) / 2);
        if (Lh - loh < n_slope_check)
          loh = std::max(hist_floor, Lh - n_slope_check);
        int nblk = Lh - loh;
        for (int k = 0; k < n_params; k++) {
          double mcb = NA_REAL, mcw = NA_REAL, est = NA_REAL;
          mc_se(k) = NA_REAL;
          if (nblk >= 2) {
            VectorXd cmean(n_chains);
            for (int c = 0; c < n_chains; c++) {
              double acc = 0;
              for (int b = loh; b < Lh; b++)
                acc += hist_chain_mean[b](c, k);
              cmean(c) = acc / (double)nblk;
            }
            est = cmean.mean();
            if (n_chains > 1) {
              double v = (cmean.array() - est).square().sum() / (n_chains - 1);
              mcb = std::sqrt(v / n_chains);
            }
            double wsum = 0;
            for (int c = 0; c < n_chains; c++) {
              double s2 = 0;
              for (int b = loh; b < Lh; b++)
                s2 += std::pow(hist_chain_mean[b](c, k) - cmean(c), 2);
              s2 /= (double)(nblk - 1);
              wsum += s2 / (double)nblk; // var of this chain's mean
            }
            mcw = std::sqrt(wsum / (double)(n_chains * n_chains));
            // Within-chain batch means, not the across-chain spread. Between-chain
            // disagreement is what R-hat already tests.
            mc_se(k) = mcw;
          }
          double info_kk = NA_REAL;
          {
            MatrixXd P = MatrixXd::Zero(n_params, n_params);
            for (int c = 0; c < n_chains; c++)
              P += opt_vec[c].get_preconditioner();
            info_kk = P(k, k) / n_chains;
          }
          conv_diag.push_back({(double)steps, (double)k, est, mcb, mcw,
                               se_stat(k), final_R_hat(k),
                               last_trend_ready ? last_tstats[k] : NA_REAL,
                               last_trend_ready ? last_rel_drift[k] : NA_REAL,
                               info_kk});
        }
      }

      if (param_conv_check && curr_batch + 1 >= n_min_batch) {
        converge =
            check_conv(win_means, win_vars, n_slope_check - 1, n_slope_check,
                       std_lim, trend_lim, par_names, print_check_info,
                       batch_steps, max_R_hat, R_hat, trend_std_conv_check,
                       R_hat_conv_check, n_chains,
                       use_std_check, trend_rel_lim, win_span_iters,
                       window_ready, curr_batch + 1, &last_rel_drift,
                       &last_conv_rhat, &last_conv_trend_std, &last_std_ratio,
                       &last_slopes, &last_tstats, &last_trend_ready);

        // Stationarity does not mean the answer is pinned down: constant-step
        // SGD settles into a distribution whose width is set by the step size.
        // Require the optimiser's own error to be small beside the statistical
        // error, se_stat = sqrt(diag(P^-1)) from the preconditioner, smoothed
        // over checkpoints. Skipped for optimisers with no preconditioner.
        std::vector<bool> mc_ok(n_params, true);
        bool mc_checked = false;
        // Precision belongs to the polish phase in stationarity mode. Leaving
        // it here would defeat the point: the search would still have to wait
        // for the Monte Carlo error of the running estimate to come down, which
        // is the very cost that averaging is supposed to absorb. se_smooth is
        // still maintained below so the polish has a precision reference.
        if (mc_se_conv_check && !stationarity_mode &&
            sgd_method == "precond_sgd") {
          for (int k = 0; k < n_params; k++) {
            if (!R_finite(se_stat(k)) || se_stat(k) <= 0)
              continue;
            se_smooth(k) = R_finite(se_smooth(k))
                               ? (0.7 * se_smooth(k) + 0.3 * se_stat(k))
                               : se_stat(k);
            if (!R_finite(mc_se(k)) || settle_left > 0)
              continue;
            mc_smooth(k) = R_finite(mc_smooth(k))
                               ? (0.7 * mc_smooth(k) + 0.3 * mc_se(k))
                               : mc_se(k);
            mc_checked = true;
            mc_ok[k] = (mc_smooth(k) <= mc_se_lim * se_smooth(k));
          }
          if (mc_checked)
            for (int k = 0; k < n_params; k++)
              converge[k] = converge[k] && mc_ok[k];
        }

        if (stationarity_mode && sgd_method == "precond_sgd") {
          for (int k = 0; k < n_params; k++) {
            if (!R_finite(se_stat(k)) || se_stat(k) <= 0) continue;
            se_smooth(k) = R_finite(se_smooth(k))
                               ? (0.7 * se_smooth(k) + 0.3 * se_stat(k))
                               : se_stat(k);
          }
        }

        // Per-parameter split-half (Geweke) stationarity: has the running
        // mean over the first half of the recent window come to agree with the
        // second half, to within its own Monte Carlo error? Scaled by the noise
        // rather than by an absolute rate, so a short window cannot pass it by
        // being uninformative.
        std::vector<bool> stationary_k(n_params, false);
        std::vector<double> stationary_z(n_params, NA_REAL);
        // Running signed / absolute half-differences, per parameter. A parameter
        // still travelling toward stationarity moves the SAME way at successive
        // checkpoints, so |sum| / sum|.| stays near 1; one that has settled but
        // wanders in a flat direction cancels and the ratio collapses.

        bool stationarity_ready = false;
        {
          const int L = (int)hist_mean.size();
          const int navail = L - hist_floor;
          if (navail >= 16) {
            stationarity_ready = true;
            // Window length from the measured correlation length. The test
            // needs enough effective samples per half: with autocorrelation
            // time tau (in sub-batches) a window of W carries only W/tau, so a
            // fixed W suits one mixing rate only. tau is the ratio of the
            // batch-means variance to the naive variance of the same mean.
            int nwin;
            if (stationarity_window > 0) {
              nwin = std::min(navail, stationarity_window);
            } else {
              const int nest = std::min(navail, 64);
              const int e0 = L - nest;
              double tau_max = 1.0;
              const int nb = std::max(2, (int)std::floor(std::sqrt((double)nest)));
              const int blen = nest / nb;
              if (blen >= 1) {
                for (int k2 = 0; k2 < n_params; k2++) {
                  double sm = 0.0, sq = 0.0;
                  for (int i2 = e0; i2 < L; ++i2) {
                    const double x = hist_mean[i2][k2]; sm += x; sq += x * x; }
                  const double mn2 = sm / (double)nest;
                  const double v_naive =
                      (sq - (double)nest * mn2 * mn2) / (double)(nest - 1) /
                      (double)nest;
                  double bs = 0.0, bq = 0.0;
                  for (int b = 0; b < nb; b++) {
                    double bm2 = 0.0;
                    for (int t = 0; t < blen; t++)
                      bm2 += hist_mean[e0 + b * blen + t][k2];
                    bm2 /= (double)blen;
                    bs += bm2; bq += bm2 * bm2;
                  }
                  const double mb2 = bs / (double)nb;
                  const double v_batch =
                      (bq - (double)nb * mb2 * mb2) / (double)(nb - 1) /
                      (double)nb;
                  if (v_naive > 0 && v_batch > 0) {
                    const double tau = v_batch / v_naive;
                    if (R_finite(tau) && tau > tau_max) tau_max = tau;
                  }
                }
              }
              if (verbose_enabled && (curr_batch % 5 == 0))
                Rcpp::Rcout << "[stat] tau_hat = " << tau_max
                            << "  navail = " << navail << std::endl;
              const double want = 2.0 * stationarity_eff * tau_max;
              nwin = (int)std::ceil(want);
              if (nwin < 16) nwin = 16;
              if (nwin > navail) stationarity_ready = false;
            }
            if (!stationarity_ready) nwin = navail;  // values unused below
            const int win0 = L - nwin;
            const int half = nwin / 2;
            const int a0 = win0, a1 = win0 + half;
            const int b0 = L - half;
            const bool have_chains =
                ((int)hist_chain_mean.size() == L) && n_chains > 1;
            for (int k = 0; k < n_params; k++) {
              // se of the half-difference. Sub-batch means are strongly
              // autocorrelated, so a naive variance understates it and the test
              // becomes impossible to pass. Take the larger of two honest
              // estimators: the spread of each chain's own half-difference
              // and batch means within the window.
              double mA = 0.0, mB = 0.0;
              for (int i2 = a0; i2 < a1; ++i2) mA += hist_mean[i2][k];
              for (int i2 = b0; i2 < L; ++i2) mB += hist_mean[i2][k];
              mA /= (double)half; mB /= (double)half;
              const double diff = std::fabs(mA - mB);

              double se_between = 0.0;
              if (have_chains) {
                double dsum = 0.0, dsq = 0.0;
                for (int c = 0; c < n_chains; c++) {
                  double cA = 0.0, cB = 0.0;
                  for (int i2 = a0; i2 < a1; ++i2) cA += hist_chain_mean[i2](c, k);
                  for (int i2 = b0; i2 < L; ++i2) cB += hist_chain_mean[i2](c, k);
                  const double dc = cA / (double)half - cB / (double)half;
                  dsum += dc; dsq += dc * dc;
                }
                const double dbar = dsum / (double)n_chains;
                const double dv =
                    (dsq - (double)n_chains * dbar * dbar) / (double)(n_chains - 1);
                if (dv > 0) se_between = std::sqrt(dv / (double)n_chains);
              }

              double se_batch = 0.0;
              {
                const int nb = std::max(2, (int)std::floor(std::sqrt((double)half)));
                const int blen = half / nb;
                if (blen >= 1 && nb >= 2) {
                  double sA = 0.0, qA = 0.0, sB = 0.0, qB = 0.0;
                  for (int b = 0; b < nb; b++) {
                    double ba = 0.0, bb = 0.0;
                    for (int t = 0; t < blen; t++) {
                      ba += hist_mean[a0 + b * blen + t][k];
                      bb += hist_mean[b0 + b * blen + t][k];
                    }
                    ba /= (double)blen; bb /= (double)blen;
                    sA += ba; qA += ba * ba; sB += bb; qB += bb * bb;
                  }
                  const double mAb = sA / nb, mBb = sB / nb;
                  const double vA = (qA - nb * mAb * mAb) / (double)(nb - 1);
                  const double vB = (qB - nb * mBb * mBb) / (double)(nb - 1);
                  if (vA > 0 && vB > 0)
                    se_batch = std::sqrt(vA / (double)nb + vB / (double)nb);
                }
              }

              const double se = std::max(se_between, se_batch);
              const bool ok = R_finite(se) && se > 0.0;

              // Magnitude, not significance. se(diff) shrinks as the chains
              // settle, so a z-test gets harder to pass the longer it runs.
              // Reference = max(se_stat, se of the half-difference). se_stat
              // alone fails at both extremes: where the data pin a parameter
              // tighter than the sampler's own noise it is unreachable.
              const double ref =
                  std::max(R_finite(se_stat(k)) ? se_stat(k) : 0.0,
                           (R_finite(se) && se > 0.0) ? se : 0.0);
              const bool have_ss = ref > 0.0;
              stationary_z[k] = have_ss ? (diff / ref) : NA_REAL;
              stationary_k[k] = have_ss && diff <= stationarity_ratio_lim * ref;
              if (stationarity_ready) {
                if (!stat_dir_init) {
                  stat_dir_hist.assign(n_params, std::deque<double>());
                  stat_flat.assign(n_params, false);
                  stat_dir_init = true;
                }
                stat_dir_hist[k].push_back(mB - mA);
                while ((int)stat_dir_hist[k].size() >
                       2 * stationarity_min_checks)
                  stat_dir_hist[k].pop_front();
              }
              // Only a parameter that is both unsettled and non-directional is
              // a flat direction. One still converging is directional, so it is
              // never flagged -- it keeps the run going, which is the point.
              if (!stationary_k[k] && stat_dir_init &&
                  (int)stat_dir_hist[k].size() >= stationarity_min_checks) {
                double dsum = 0.0, dabs = 0.0;
                for (double v : stat_dir_hist[k]) {
                  dsum += v; dabs += std::fabs(v); }
                if (dabs <= 0.0) continue;
                const double dir = std::fabs(dsum) / dabs;
                if (dir < stationarity_dir_lim && !stat_flat[k]) {
                  stat_flat[k] = true;
                  if (verbose_enabled) {
                    ngme_io::err() << "[stationarity] "
                                << (k < (int)par_names.size() ? par_names[k]
                                                             : std::string("?"))
                                << ": flat direction (moves "
                                << diff / se_stat(k)
                                << " x se_stat but cancels, |sum|/sum|.| = "
                                << dir << "); excluded from the stopping rule."
                                << std::endl;
                    R_FlushConsole();
                  }
                }
              }
            }
          }
        }
        if (stationarity_mode) {
          for (int k = 0; k < n_params; k++)
            if (!(stationarity_ready &&
                  (stationary_k[k] || (stat_dir_init && stat_flat[k]))))
              converge[k] = false;
        }
        // Compact progress for long runs. Printed AFTER every clause has been
        // applied to converge[], and reporting the statistic that decides.
        if (verbose_enabled && !print_check_info) {
          int n_ok = 0, i_rhat = 0, i_z = 0;
          for (int i = 0; i < n_params; i++) {
            if (converge[i]) n_ok++;
            if (final_R_hat(i) > final_R_hat(i_rhat)) i_rhat = i;
            if (R_finite(stationary_z[i]) &&
                (!R_finite(stationary_z[i_z]) ||
                 stationary_z[i] > stationary_z[i_z]))
              i_z = i;
          }
          auto nm = [&](int i) {
            return (i < (int)par_names.size()) ? par_names[i] : std::string("?");
          };
          ngme_io::err() << "[iter " << steps << "] " << n_ok << "/" << n_params
                      << " converged";
          if (!stationarity_mode) {
            // Forecast only where precision gates the search. Under
            // stationarity it gates the polish instead, and mc_smooth is not
            // maintained here, so there is nothing to forecast from.
            double worst_ratio = 0.0;
            for (int k = 0; k < n_params; k++)
              if (R_finite(mc_smooth(k)) && R_finite(se_smooth(k)) &&
                  se_smooth(k) > 0)
                worst_ratio = std::max(
                    worst_ratio, mc_smooth(k) / (mc_se_lim * se_smooth(k)));
            if (worst_ratio > 1.0)  // mc_se ~ 1/sqrt(n), so this is a forecast
              ngme_io::err() << "  [precision needs ~"
                          << (long)std::ceil(steps * worst_ratio * worst_ratio)
                          << " iters]";
          }
          if (n_ok < n_params) {
            ngme_io::err() << "  worst: R_hat " << std::fixed
                        << std::setprecision(3) << final_R_hat(i_rhat) << " ("
                        << nm(i_rhat) << ")";
            if (stationarity_mode && R_finite(stationary_z[i_z]))

              ngme_io::err() << ", move/ref " << std::setprecision(2)
                          << stationary_z[i_z] << " (" << nm(i_z) << ", pass <= "
                          << std::setprecision(2) << stationarity_ratio_lim
                          << ")";
            else if (last_trend_ready)
              ngme_io::err() << ", drift/100 " << std::setprecision(2)
                          << (100.0 * last_rel_drift[0]) << "%";
          }
          ngme_io::err() << std::endl;
          // Current parameter values reported on the optimiser's theta scale.
          ngme_io::err() << "           theta:";
          for (int i = 0; i < n_params; i++)
            ngme_io::err() << " " << nm(i) << "=" << std::fixed
                           << std::setprecision(4) << means(curr_batch, i);
          ngme_io::err() << std::defaultfloat << std::setprecision(6)
                         << std::endl;
          R_FlushConsole();
        }


        if (schedule_auto_start && !schedule_armed) {
          const int L = (int)hist_mean.size();
          const int navail = L - hist_floor;
          if (navail >= 8) {
            const int half = navail / 2;
            const int a0 = hist_floor, a1 = hist_floor + half;
            const int b0 = L - half;
            bool stationary = true;
            for (int k = 0; k < n_params && stationary; k++) {
              double mA = 0.0, mB = 0.0;
              for (int i = a0; i < a1; ++i) mA += hist_mean[i][k];
              for (int i = b0; i < L; ++i) mB += hist_mean[i][k];
              mA /= (double)half;
              mB /= (double)half;
              double vA = 0.0, vB = 0.0;
              for (int i = a0; i < a1; ++i) {
                const double d = hist_mean[i][k] - mA; vA += d * d;
              }
              for (int i = b0; i < L; ++i) {
                const double d = hist_mean[i][k] - mB; vB += d * d;
              }
              // batch-means standard error of each half-average
              vA /= (double)half * (double)(half - 1);
              vB /= (double)half * (double)(half - 1);
              const double se = std::sqrt(vA + vB);
              if (!R_finite(se) || !(se > 0.0) ||
                  std::fabs(mA - mB) > schedule_arm_z * se)
                stationary = false;
            }
            if (stationary)
              schedule_ready_count++;
            else
              schedule_ready_count = 0;
            if (schedule_ready_count >= n_conv_batch) {
              schedule_armed = true;
              for (int c = 0; c < n_chains; c++)
                opt_vec[c].set_schedule_burnin(steps);
              if (verbose_enabled)
                Rcpp::Rcout << "[schedule] Polyak average stationary at "
                            << "iteration " << steps
                            << "; step-size decay starts here" << std::endl;
            }
          }
        }

        bool all_params_ok =
            std::find(begin(converge), end(converge), false) == end(converge);

        // Stationary but imprecise: shrink eta rather than iterate. mc_se
        // falls as 1/sqrt(n) around a distribution whose width is set by eta,
        // so narrowing it (width ~ sqrt(eta)) is the standard remedy. Only
        // after the transient is over, and capped so a hopeless fit ends.
        if (settle_left > 0)
          settle_left--;
        if (!all_params_ok && mc_checked && settle_left == 0 &&
            n_stepsize_decays < max_stepsize_decays) {
          bool stationary = true;
          for (int k = 0; k < n_params; k++) {
            if (R_hat_conv_check && !last_conv_rhat[k]) stationary = false;
            if (trend_std_conv_check &&
                (!last_trend_ready || !last_conv_trend_std[k])) stationary = false;
          }
          if (stationary) {
            n_stepsize_decays++;
            stepsize_decay_scale *= stepsize_decay_precision_gamma;
            // everything measured so far belongs to the old step size
            hist_floor = (int)hist_mean.size();
            settle_left = n_settle_checks;
            mc_smooth = VectorXd::Constant(n_params, NA_REAL);
            for (int c = 0; c < n_chains; c++) {
              opt_vec[c].set_stepsize_decay_enabled(true);
              opt_vec[c].set_stepsize_decay_scale(stepsize_decay_scale);
            }
            consec_converged = 0;
            if (verbose_enabled) {
              double worst = 0.0;
              int kw = 0;
              for (int k = 0; k < n_params; k++)
                if (R_finite(mc_se(k)) && R_finite(se_smooth(k)) &&
                    se_smooth(k) > 0 && mc_se(k) / se_smooth(k) > worst) {
                  worst = mc_se(k) / se_smooth(k);
                  kw = k;
                }
              ngme_io::err() << "[iter " << steps
                          << "] stationary but imprecise (mc_se/se = "
                          << std::fixed << std::setprecision(2) << worst << " on "
                          << (kw < (int)par_names.size() ? par_names[kw] : "?")
                          << "); stepsize x " << stepsize_decay_precision_gamma
                          << " -> scale " << stepsize_decay_scale << std::endl;
            }
          }
        }
        // Require the criteria to hold on n_conv_batch consecutive checkpoints.
        // A single passing checkpoint is weak evidence: R-hat bounces around
        // its threshold from batch to batch even long after the chains mix.
        consec_converged = all_params_ok ? consec_converged + 1 : 0;
        batch_converged = (consec_converged >= n_conv_batch);
        if (batch_converged)
          converged_by_param = true;

        // Trend-triggered step-size reduction.
        if (stepsize_decay_on_trend && last_trend_ready &&
            trend_decay_cooldown <= 0) {
          bool no_drift = true;
          for (int i = 0; i < n_params; i++)
            if (!(R_finite(last_rel_drift[i]) &&
                  last_rel_drift[i] <= trend_rel_lim))
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
  // one passing check is a coin flip on a noisy statistic; require two
  int polish_precise_streak = 0;
  // steps grows during the polish phase too, so remember where the stopping
  // rule actually fired
  const int steps_at_convergence = steps;
  // Only polish a converged run. The "no longer than the search" cap applies
  // to drift mode only: under stationarity the search is deliberately short and
  // the averaging is what buys precision, so that cap would bound the precision
  // by the length of the search.
  const int polish_budget =
      all_converge ? (stationarity_mode ? polish_iterations
                                        : std::min(polish_iterations, steps))
                   : 0;
  // How often to ask whether the average is precise enough to stop.
  const int polish_check_every = std::max(20, batch_steps);
  // Per-chain batch sums over the polish window, so the Monte Carlo error of
  // the average can be estimated from batch means as well as from the spread
  // between chains.
  const int polish_batch_len = std::max(10, polish_check_every / 4);
  // Step floor from the chains' own spread.  The polish decays the step to
  // remove the O(eta) bias, but averaging only reduces error while the chains
  // still mix.
  MatrixXd polish_cur = MatrixXd::Zero(n_chains, n_params);
  MatrixXd polish_prev = MatrixXd::Zero(n_chains, n_params);
  MatrixXd polish_incr_sq = MatrixXd::Zero(n_chains, n_params);
  int polish_incr_n = 0;
  bool polish_prev_set = false;
  double polish_floor_scale = 0.0;  // 0 = no floor in force
  std::vector<bool> precision_limited(n_params, false);
  std::vector<std::vector<Eigen::VectorXd>> polish_batch_means(n_chains);
  MatrixXd polish_batch_sum = MatrixXd::Zero(n_chains, n_params);
  int polish_batch_fill = 0;
  if (polish_budget > 0 && n_chains > 0) {
    double polish_scale = stepsize_decay_scale * polish_stepsize_factor;
    double polish_cur_scale = polish_scale;
    for (i = 0; i < n_chains; i++) {
      opt_vec[i].set_stepsize_decay_enabled(true);
      opt_vec[i].set_stepsize_decay_scale(polish_scale);
    }
    std::atomic<bool> polish_failed(false);
    std::string polish_error;
    for (int step = 0; step < polish_budget; step++) {
      // Normalised so the decay starts at exactly the polish scale, matching
      // the schedule used in the search phase. No floor here: the whole point
      // is to let eta go to zero so the O(eta) bias goes with it, and the
      // averaging is what controls the variance.
      if (polish_decay_alpha > 0.0) {
        double sc =
            polish_scale *
            std::pow(1.0 + (double)step / polish_decay_t0, -polish_decay_alpha);
        // Never below the spread-derived floor, and the floor itself is capped
        // at the step the polish started from, so the worst case is that the
        // polish degrades to constant-step averaging and never an increase.
        if (polish_floor_scale > 0.0 && sc < polish_floor_scale)
          sc = polish_floor_scale;
        for (int c = 0; c < n_chains; c++)
          opt_vec[c].set_stepsize_decay_scale(sc);
        polish_cur_scale = sc;
      }
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
          {
            polish_sum.row(i) += param;
            polish_batch_sum.row(i) += param;
            polish_cur.row(i) = param.transpose();
          }
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
      if (polish_prev_set) {
        polish_incr_sq.array() += (polish_cur - polish_prev).array().square();
        polish_incr_n++;
      }
      polish_prev = polish_cur;
      polish_prev_set = true;
      if (++polish_batch_fill >= polish_batch_len) {
        for (int c = 0; c < n_chains; c++)
          polish_batch_means[c].push_back(
              (polish_batch_sum.row(c) / (double)polish_batch_fill).transpose());
        polish_batch_sum.setZero();
        polish_batch_fill = 0;
      }

      // Stop when the average is precise. The chains are independent, so the
      // Monte Carlo error of the reported estimate is the spread of
      // those per-chain averages divided by sqrt(n_chains). Polishing until
      // that sits under mc_se_lim * se_stat asks for the Monte Carlo error to
      // be a fraction of the statistical error, which is the point at which
      // further iterations stop buying anything the data can support.
      if (n_chains > 1 && polish_done >= polish_check_every &&
          (polish_done % polish_check_every == 0)) {
        // Classify each parameter by whether averaging can still work.
        // r_k = delta_k * sqrt(H) / S_k, dimensionless: can this chain's random
        // walk cross the current between-chain spread within the horizon?
        // delta_k is measured after the preconditioner, so no per-parameter
        // scale knowledge is needed and the same rule applies to every model.
        if (polish_incr_n > 0) {
          const double H = (double)std::max(1, polish_budget);
          for (int k = 0; k < n_params; k++) {
            if (!R_finite(se_stat(k)) || se_stat(k) <= 0) continue;
            // delta: rms per-iteration increment, pooled over chains. The
            // increment of a random walk has variance 2*sigma^2, hence the /2.
            double dsq = 0.0;
            for (int c = 0; c < n_chains; c++) dsq += polish_incr_sq(c, k);
            const double delta =
                std::sqrt(dsq / (double)(n_chains * polish_incr_n) / 2.0);
            // S: spread of the chains' current positions
            double m = 0.0;
            for (int c = 0; c < n_chains; c++) m += polish_cur(c, k);
            m /= (double)n_chains;
            double sv = 0.0;
            for (int c = 0; c < n_chains; c++) {
              const double d = polish_cur(c, k) - m; sv += d * d; }
            const double S = (n_chains > 1)
                                 ? std::sqrt(sv / (double)(n_chains - 1)) : 0.0;
            if (!(S > 0.0) || !R_finite(delta) || delta <= 0.0) continue;
            const double r = delta * std::sqrt(H) / S;
            // delta scales linearly with the step, so r at the polish's
            // STARTING step is r rescaled by that ratio.
            const double r_at_start =
                (polish_cur_scale > 0.0) ? r * (polish_scale / polish_cur_scale)
                                         : r;
            if (r_at_start < 1.0) {
              // REGIME 3: unreachable even without any decay. No schedule fixes
              // this -- the chains cannot cross their own spread at full step.
              // Report it instead of letting it burn the whole budget.
              if (!precision_limited[k] && verbose_enabled) {
                ngme_io::err() << "[polish] "
                            << (k < (int)par_names.size() ? par_names[k] : std::string("?"))
                            << ": precision-limited (chains cannot cross their"
                            << " spread at full step, r = " << r_at_start
                            << "); excluded from the stopping rule."
                            << std::endl;
                R_FlushConsole();
              }
              precision_limited[k] = true;
            }
          }
          // Constant lower bound, a fraction of the polish-start step. No
          // bound and the chains freeze; a bound that feeds back on r and the
          // step never decays. A constant does neither.
          if (polish_floor_scale <= 0.0 && polish_floor_frac > 0.0) {
            polish_floor_scale = polish_floor_frac * polish_scale;
            if (verbose_enabled) {
              Rcpp::Rcout << "[polish] step bounded below at "
                          << polish_floor_scale << " (" << polish_floor_frac
                          << " x the polish-start step)" << std::endl;
              R_FlushConsole();
            }
          }
          polish_incr_sq.setZero();
          polish_incr_n = 0;
        }
        // Exclusions are for the exceptional parameter. If most are excluded
        // the rule tests almost nothing.
        int n_checkable = 0, n_excluded = 0;
        for (int k = 0; k < n_params; k++) {
          if (!R_finite(se_stat(k)) || se_stat(k) <= 0) continue;
          if (precision_limited[k]) n_excluded++; else n_checkable++;
        }
        bool precise = (n_checkable > n_excluded);
        if (!precise && verbose_enabled && polish_precise_streak == 0)
          ngme_io::err() << "[polish] " << n_excluded << " of "
                      << (n_checkable + n_excluded)
                      << " parameters are precision-limited: too many to"
                      << " certify precision, polishing to the cap."
                      << std::endl;
        for (int k = 0; k < n_params && precise; k++) {
          if (!R_finite(se_stat(k)) || se_stat(k) <= 0) continue;
          if (precision_limited[k]) continue;  // regime 3: cannot be met
          double m = 0.0;
          for (int c = 0; c < n_chains; c++)
            m += polish_sum(c, k) / (double)polish_done;
          m /= (double)n_chains;
          double v = 0.0;
          for (int c = 0; c < n_chains; c++) {
            const double d = polish_sum(c, k) / (double)polish_done - m;
            v += d * d;
          }
          v /= (double)(n_chains - 1);
          double mc_se = std::sqrt(v / (double)n_chains);

          // Batch-means estimate of the same quantity, pooled over chains.
          // var(chain average) ~ var(batch means)/B, and the estimate reported
          // is the mean over chains, so its variance is that over n_chains.
          long df = 0; double pooled = 0.0;
          for (int c = 0; c < n_chains; c++) {
            const int B = (int)polish_batch_means[c].size();
            if (B < 2) continue;
            double bm = 0.0;
            for (int b = 0; b < B; b++) bm += polish_batch_means[c][b](k);
            bm /= (double)B;
            double bv = 0.0;
            for (int b = 0; b < B; b++) {
              const double d = polish_batch_means[c][b](k) - bm; bv += d * d; }
            pooled += bv / (double)(B - 1) / (double)B;  // var of this chain's mean
            df += (B - 1);
          }
          if (df >= 4) {
            const double mc_se_bm =
                std::sqrt(pooled / (double)(n_chains * n_chains));
            // Take the larger of the two: they estimate the same quantity, and
            // being wrong in the conservative direction only costs iterations,
            // whereas being wrong the other way reports a precision we do not
            // have.
            mc_se = std::max(mc_se, mc_se_bm);
          }
          if (!(mc_se <= mc_se_lim * se_stat(k))) precise = false;

          // Has the reported average stopped moving? A small mc_se only says
          // the chains agree about it. The relative target is most lenient
          // where a parameter is worst identified.
          if (precise) {
            double dnum = 0.0, dvar = 0.0; int dchains = 0;
            for (int c = 0; c < n_chains; c++) {
              const int B = (int)polish_batch_means[c].size();
              if (B < 4) continue;
              const int h = B / 2;
              double m1 = 0.0, m2 = 0.0;
              for (int b = 0; b < h; b++) m1 += polish_batch_means[c][b](k);
              for (int b = B - h; b < B; b++) m2 += polish_batch_means[c][b](k);
              m1 /= (double)h; m2 /= (double)h;
              double v = 0.0, bm = 0.0;
              for (int b = 0; b < B; b++) bm += polish_batch_means[c][b](k);
              bm /= (double)B;
              for (int b = 0; b < B; b++) {
                const double d = polish_batch_means[c][b](k) - bm; v += d * d; }
              v /= (double)(B - 1);
              dnum += (m2 - m1);
              dvar += 2.0 * v / (double)h;  // var of the difference of halves
              dchains++;
            }
            if (dchains > 0) {
              const double diff = std::abs(dnum / (double)dchains);
              const double se_d = std::sqrt(dvar) / (double)dchains;
              const double ref_d =
                  std::max(R_finite(se_stat(k)) ? se_stat(k) : 0.0,
                           (R_finite(se_d) && se_d > 0.0) ? se_d : 0.0);
              if (ref_d > 0.0 && diff > mc_se_lim * ref_d)
                precise = false;  // the average is still moving
            }
          }
        }
        polish_precise_streak = precise ? (polish_precise_streak + 1) : 0;
        if (precise && polish_precise_streak >= 2) {
          if (verbose_enabled)
            ngme_io::err() << "Polish reached target precision after "
                        << polish_done << " iterations" << std::endl;
          break;
        }
      }
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
  {
    // conv_diag -> a numeric matrix R can read straight into a data.frame
    const int ncol = 10;
    Rcpp::NumericMatrix cd((int)conv_diag.size(), ncol);
    for (size_t r = 0; r < conv_diag.size(); r++)
      for (int c = 0; c < ncol; c++)
        cd(r, c) = conv_diag[r][c];
    Rcpp::colnames(cd) = Rcpp::CharacterVector::create(
        "iteration", "param", "estimate", "mc_se_between", "mc_se_batch",
        "se_stat", "R_hat", "t_stat", "drift_per100", "info_kk");
    outputs.attr("conv_diag") = cd;
  }
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
           bool use_std_check, double trend_rel_lim,
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

  // Computed whenever the window is ready. These are the drift diagnostics reported in
  // attr(fit, "conv_diag"); leaving them uncomputed under a different criterion
  // meant conv_diag reported drift_per100 = 0 and t_stat = 0, which reads as
  // "no drift" when it means "not measured".
  if (trend_ready) {
    for (int i = 0; i < n_params; i++) {
      std_ratio[i] = std::sqrt(vars(curr_batch, i)) /
                     (std::abs(means(curr_batch, i)) + 1e-5);
      // The raw coefficient of variation is not a convergence statistic: it
      // conflates "the chains disagree" with "this parameter is imprecisely
      // determined", so a weakly identified parameter can never pass no
      // matter how well mixed the chains are. Off by default as R-hat already
      // measures between- vs within-chain disagreement in a scale-free way.
      if (trend_std_conv_check && use_std_check && std_ratio[i] > std_lim) {
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

      // Drift as a rate per 100 iterations, relative to the parameter's scale.
      double span_pts = (double)(n_slope_check - 1);
      double scale = std::abs(means(curr_batch, i)) +
                     std::sqrt(std::max(0.0, vars(curr_batch, i))) + 1e-10;
      double per_iter = (win_span_iters > 0.0)
                            ? (std::abs(beta(1)) * span_pts / win_span_iters)
                            : 0.0;
      rel_drift[i] = 100.0 * per_iter / scale; // fraction of scale per 100 iters

      // Drift-rate test only.
      if (trend_std_conv_check && rel_drift[i] > trend_rel_lim) {
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
    // curr_batch indexes the regression window, so the checkpoint number for the reader
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
