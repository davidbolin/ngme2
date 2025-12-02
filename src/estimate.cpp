#ifdef USEMKL
#define EIGEN_USE_MKL_ALL
#endif

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include "include/timer.h"
#include "ngme.h"
#include "optimizer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
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
#ifdef _OPENMP
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

  const bool verbose_enabled = control_opt["verbose"];

  VectorXi num_threads = Rcpp::as<VectorXi>(control_opt["num_threads"]);
  int n_threads_chain = num_threads[0];

  omp_set_max_active_levels(2);
  omp_set_num_threads(num_threads[0] * num_threads[1]);

  // init model and optimizer
  vector<std::shared_ptr<Ngme>> ngmes;
  vector<Ngme_optimizer> opt_vec;
  int i = 0;
  for (i = 0; i < n_chains; i++) {
    double sd = (i == 0) ? 0 : start_sd;

    // Not thread-safe using Rcpp::List to init optimizer
    ngmes.push_back(std::make_shared<Ngme>(R_ngme, seed + i, sampling_strategy,
                                           num_threads[1], sd));
    opt_vec.push_back(Ngme_optimizer(control_opt, ngmes[i]));
    opt_vec.back().set_pflug_conv_check(pflug_conv_check);
    if (verbose_enabled && i > 0) {
      opt_vec.back().set_verbose(false);
    }
  }

  std::vector<std::string> par_names = ngmes[0]->get_par_names();

// burn in period
#pragma omp parallel for schedule(static) num_threads(n_threads_chain)
  for (i = 0; i < n_chains; i++)
    (ngmes[i])->burn_in(burnin);

  int n_params = ngmes[0]->get_n_params();
  MatrixXd means(n_batch, n_params);
  MatrixXd vars(n_batch, n_params);

  // for Gelman-Rubin statistic
  MatrixXd batch_sum(n_chains, n_params);
  MatrixXd batch_sq_sum(n_chains, n_params);
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

  int curr_batch = 0;

  while (steps < iterations && !all_converge) {
    MatrixXd mat(n_chains, n_params);
    batch_sum.setZero();
    batch_sq_sum.setZero();

    // Store max_pflug_sum before batch
    std::vector<double> max_pflug_sum_before(n_chains);
    if (pflug_conv_check) {
      for (int i = 0; i < n_chains; i++) {
        max_pflug_sum_before[i] = opt_vec[i].get_max_pflug_sum();
      }
    }

    // Run batch_steps iterations using unified SGD step (computes grad and,
    // optionally, precond)
    for (int step = 0; step < batch_steps; step++) {
#pragma omp parallel for schedule(static) num_threads(n_threads_chain)
      for (i = 0; i < n_chains; i++) {
        // Compute one SGD step; decide whether to compute preconditioner this
        // iter We compute it every iter if compute_precond_each_iter, else let
        // the optimizer refresh at needed cadence
        VectorXd param = opt_vec[i].sgd(0.1, // eps
                                        1,   // one step per loop
                                        max_relative_step, max_absolute_step,
                                        compute_precond_each_iter);
        // if (step == batch_steps - 1) {
#pragma omp critical
        {
          mat.row(i) = param;
          batch_sum.row(i) += param;
          batch_sq_sum.row(i) += param.array().square().matrix();
        }
        // }
      }
    }
    steps += batch_steps;

    // compute mean and variance
    means.row(curr_batch) = mat.colwise().mean();
    for (int k = 0; k < n_params; k++) {
      vars(curr_batch, k) =
          (mat.col(k).array() - means(curr_batch, k)).square().sum() /
          (n_chains - 1);
    }

    if (n_chains > 1) {
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

      Eigen::RowVectorXd var_hat = ((double)(n - 1) / n) * W + (1.0 / n) * B;
      R_hat = (var_hat.array() / W.array()).sqrt().transpose();
      final_R_hat = R_hat;

      if (curr_batch + 1 >= n_min_batch) {
        converge =
            check_conv(means, vars, curr_batch, n_slope_check, std_lim,
                       trend_lim, par_names, print_check_info, batch_steps,
                       max_R_hat, R_hat, trend_std_conv_check, R_hat_conv_check,
                       &last_conv_rhat, &last_conv_trend_std, &last_std_ratio,
                       &last_slopes, &last_trend_ready);
        all_converge =
            std::find(begin(converge), end(converge), false) == end(converge);
        if (all_converge)
          converged_by_param = true;
      } else {
        all_converge = false;
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
          all_converge = true;
          pflug_triggered = true;
          if (verbose_enabled)
            Rcpp::Rcout << "Pflug diagnostic satisfied: pflug_sum < "
                        << pflug_alpha << " * max_pflug_sum for all chains.\n"
                        << std::endl;
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

#else // No parallel chain
  Ngme ngme(R_ngme, seed, sampling_strategy);
  Ngme_optimizer opt(control_opt, std::make_shared<Ngme>(ngme));
  opt.set_pflug_conv_check(pflug_conv_check);
  opt.sgd(0.1, iterations, max_relative_step, max_absolute_step,
          compute_precond_each_iter);
  // estimation done, posterior sampling
  // ngme.sampling(10, true);
  outputs.push_back(ngme.output());
  if (store_traj) {
    opt.record_current_state();
    trajs_chains.push_back(opt.get_trajs());
  }
#endif

  std::cout << "Estimation ends." << std::endl;
  std::cout << "Total time of the estimation is (s): "
            << since(timer).count() / 1000 << std::endl;

  if (store_traj) {
    outputs.attr("opt_traj") = trajs_chains;
  } else {
    outputs.attr("opt_traj") = R_NilValue;
  }
#ifdef _OPENMP
  if (n_chains > 1)
    outputs.attr("R_hat") = final_R_hat;
#endif
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
