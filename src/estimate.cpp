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
                             std::vector<bool> *conv_rhat_out,
                             std::vector<bool> *conv_trend_std_out,
                             std::vector<double> *std_ratio_out,
                             std::vector<double> *slopes_out,
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
  const bool pflug_conv_check =
      control_opt.containsElementNamed("pflug_conv_check")
          ? Rcpp::as<bool>(control_opt["pflug_conv_check"])
          : false;
  const double pflug_alpha = control_opt.containsElementNamed("pflug_alpha")
                                 ? Rcpp::as<double>(control_opt["pflug_alpha"])
                                 : 1.0;

  Rcpp::List output = R_NilValue;

  vector<vector<VectorXd>> trajs_chains;

  auto timer = std::chrono::steady_clock::now();

  Rcpp::List outputs;
  const bool verbose_enabled = control_opt["verbose"];
  const std::string sgd_method =
      Rcpp::as<std::string>(control_opt["sgd_method"]);

  bool stepsize_decay_enabled = false;
  int stepsize_decay_patience = 0;
  double stepsize_decay_gamma = 1.0;
  double stepsize_decay_min_delta = 0.0;
  int stepsize_decay_warmup = 0;
  if (control_opt.containsElementNamed("stepsize_decay")) {
    std::string decay = Rcpp::as<std::string>(control_opt["stepsize_decay"]);
    stepsize_decay_enabled = (decay == "grad_norm_plateau");
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

  // init model and optimizer
  vector<std::shared_ptr<Ngme>> ngmes;
  vector<Ngme_optimizer> opt_vec;
  int i = 0;
  for (i = 0; i < n_chains; i++) {
    double sd = (i == 0) ? 0 : start_sd;

    // Not thread-safe using Rcpp::List to init optimizer
    ngmes.push_back(std::make_shared<Ngme>(R_ngme, seed + i, sampling_strategy,
                                           num_threads[1], sd));
    opt_vec.push_back(Ngme_optimizer(control_opt, ngmes[i], seed + i));
    opt_vec.back().set_pflug_conv_check(pflug_conv_check);
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
  // keep last-diagnostics for reporting
  std::vector<bool> last_conv_rhat(n_params, false);
  std::vector<bool> last_conv_trend_std(n_params, false);
  std::vector<double> last_std_ratio(n_params, 0.0);
  std::vector<double> last_slopes(n_params, 0.0);
  bool last_trend_ready = false;
  bool pflug_triggered = false;
  bool converged_by_param = false;
  std::vector<double> last_pflug_sum(n_chains, 0.0);
  std::vector<double> last_pflug_max(n_chains, 0.0);

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

  int curr_batch = 0;
  double stepsize_decay_scale = 1.0;
  double stepsize_decay_prev_norm = std::numeric_limits<double>::infinity();
  int stepsize_decay_bad_epochs = 0;

  // Enable convergence only when at least one diagnostic is requested.
  const bool param_conv_check = trend_std_conv_check || R_hat_conv_check;
  const bool any_conv_check = param_conv_check || pflug_conv_check;

  while (steps < iterations && !all_converge) {
    MatrixXd mat(n_chains, n_params);
    batch_sum.setZero();
    batch_sq_sum.setZero();
    for (int b = 0; b < n_sub_batch; ++b) {
      sub_batch_sum[b].setZero();
      sub_batch_count[b] = 0;
    }

    // Store max_pflug_sum before batch
    std::vector<double> max_pflug_sum_before(n_chains);
    if (pflug_conv_check) {
      for (int i = 0; i < n_chains; i++) {
        max_pflug_sum_before[i] = opt_vec[i].get_max_pflug_sum();
      }
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

      if (param_conv_check && curr_batch + 1 >= n_min_batch) {
        converge =
            check_conv(means, vars, curr_batch, n_slope_check, std_lim,
                       trend_lim, par_names, print_check_info, batch_steps,
                       max_R_hat, R_hat, trend_std_conv_check, R_hat_conv_check,
                       &last_conv_rhat, &last_conv_trend_std, &last_std_ratio,
                       &last_slopes, &last_trend_ready);
        batch_converged =
            std::find(begin(converge), end(converge), false) == end(converge);
        if (batch_converged)
          converged_by_param = true;
      }

      // 2. if some parameter converge, stop compute gradient, or slow down the
      // gradient. if (auto_stop)
      //     for (int i=0; i < n_chains; i++) {
      //         ngmes[i]->check_converge(converge);
      //     }
      // Pflug diagnostic check for parallel chains
      if (pflug_conv_check && curr_batch + 1 >= n_min_batch) {
        bool pflug_converged = true;
        for (int i = 0; i < n_chains; i++) {
          last_pflug_sum[i] = opt_vec[i].get_pflug_sum();
          last_pflug_max[i] = opt_vec[i].get_max_pflug_sum();
        }
        for (int i = 0; i < n_chains; i++) {
          double curr_sum = opt_vec[i].get_pflug_sum();
          double curr_max = opt_vec[i].get_max_pflug_sum();
          double threshold = pflug_alpha * curr_max;
          // Require curr_sum < threshold; if curr_max == 0, treat as not yet
          // converged
          if (curr_max <= 0 || curr_sum >= threshold) {
            pflug_converged = false;
            break;
          }
        }
        if (pflug_converged) {
          batch_converged = true;
          pflug_triggered = true;
          if (verbose_enabled)
            Rcpp::Rcout << "Pflug diagnostic satisfied: pflug_sum < "
                        << pflug_alpha << " * max_pflug_sum for all chains.\n"
                        << std::endl;
        }
      }
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

  // generate outputs
  for (i = 0; i < n_chains; i++) {
    // use last batch average parameter instead of last step value
    Eigen::VectorXd avg_param = (batch_sum.row(i) / batch_steps).transpose();
    ngmes[i]->set_parameter_and_update(avg_param, false);

    outputs.push_back(ngmes[i]->output());
    if (store_traj) {
      opt_vec[i].record_current_state();
      trajs_chains.push_back(opt_vec[i].get_trajs());
    }
  }
  if (all_converge)
    std::cout << "Reach convergence in " << steps << " iterations."
              << std::endl;

  if (all_converge && n_chains > 1) {
    bool converged_by_pflug_only = pflug_triggered && !converged_by_param;
    std::cout << "Convergence criteria summary:\n";

    if (converged_by_pflug_only) {
      std::cout << "  - Pflug diagnostic satisfied (pflug_sum < " << pflug_alpha
                << " * max_pflug_sum for all chains)\n";
      std::cout << "  Per-chain Pflug stats (sum / max):\n";
      for (int i = 0; i < n_chains; i++) {
        std::cout << "    * chain " << i + 1 << ": " << std::fixed
                  << std::setprecision(4) << last_pflug_sum[i] << " / "
                  << last_pflug_max[i] << "\n";
      }
    } else {
      if (R_hat_conv_check) {
        std::cout << "  - R_hat threshold: max_R_hat = " << max_R_hat << "\n";
      }
      if (trend_std_conv_check) {
        std::cout << "  - Trend/Std thresholds: std_lim = " << std_lim
                  << ", trend_lim = " << trend_lim
                  << ", window (stop points) = " << n_slope_check << "\n";
      }

      std::cout << "  Per-parameter status (R_hat | std/mean | slope):\n";
      for (int i = 0; i < n_params; i++) {
        std::cout << "    * " << par_names[i] << ": ";
        bool printed_any = false;
        if (R_hat_conv_check) {
          std::cout << "R_hat=" << std::fixed << std::setprecision(3)
                    << final_R_hat(i)
                    << (last_conv_rhat[i] ? " (ok)" : " (fail)");
          printed_any = true;
        }
        if (trend_std_conv_check && last_trend_ready) {
          if (printed_any)
            std::cout << "; ";
          std::cout << "std/mean=" << std::setprecision(3) << last_std_ratio[i]
                    << (last_conv_trend_std[i] ? " (ok)" : " (fail)") << ", ";
          std::cout << "slope=" << std::setprecision(3) << last_slopes[i]
                    << (std::abs(last_slopes[i]) <= trend_lim ? " (ok)"
                                                              : " (fail)");
        }
        std::cout << "\n";
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
           bool trend_std_conv_check, bool R_hat_conv_check,
           std::vector<bool> *conv_rhat_out,
           std::vector<bool> *conv_trend_std_out,
           std::vector<double> *std_ratio_out, std::vector<double> *slopes_out,
           bool *trend_ready_out) {
  int n_params = means.cols();
  std::vector<bool> conv(n_params, false);

  std::vector<bool> conv_rhat(n_params, true);
  std::vector<bool> conv_trend_std(n_params, true);
  std::vector<double> std_ratio(n_params, 0.0);
  std::vector<double> slopes(n_params, 0.0);

  bool trend_ready = (curr_batch + 1 >= n_slope_check);

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
      if (std_ratio[i] > std_lim) {
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
      if (std::abs(beta(1)) > trend_lim) {
        conv_trend_std[i] = false;
      }
    }
  }

  for (int i = 0; i < n_params; i++) {
    bool passed = false;
    if (trend_std_conv_check && trend_ready && conv_trend_std[i])
      passed = true;
    if (R_hat_conv_check && conv_rhat[i])
      passed = true;
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
  if (trend_ready_out)
    *trend_ready_out = trend_ready;

  if (print_check_info) {
    std::cout << "\nstop " << curr_batch + 1 << ":\n";

    const int label_width = 11; // width for the row label (e.g., "R_hat:")
    const int col_width = 10;   // width for each value/parameter name
    int line_width = label_width + n_params * (col_width + 1);

    auto print_separator = [&]() {
      std::cout << std::string(line_width, '-') << "\n";
    };

    print_separator();

    std::cout << std::setw(label_width) << std::left << "Param:";
    for (const auto &name : par_names) {
      std::cout << " " << std::setw(col_width) << std::left << name;
    }
    std::cout << "\n";

    print_separator();

    if (R_hat_conv_check) {
      std::cout << std::setw(label_width) << std::left << "R_hat:";
      for (int i = 0; i < n_params; i++) {
        std::cout << " " << std::setw(col_width) << std::fixed
                  << std::setprecision(3) << std::left << R_hat(i);
      }
      std::cout << "\n";
    }

    if (trend_std_conv_check && trend_ready) {
      std::cout << std::setw(label_width) << std::left << "std/mean:";
      for (int i = 0; i < n_params; i++) {
        std::cout << " " << std::setw(col_width) << std::fixed
                  << std::setprecision(3) << std::left << std_ratio[i];
      }
      std::cout << "\n";

      std::cout << std::setw(label_width) << std::left << "trend:";
      for (int i = 0; i < n_params; i++) {
        std::cout << " " << std::setw(col_width) << std::fixed
                  << std::setprecision(3) << std::left << slopes[i];
      }
      std::cout << "\n";
    }

    print_separator();
    std::cout << "\n";
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
