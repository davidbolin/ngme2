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
                             bool, int, double, const VectorXd &);

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

  Rcpp::List output = R_NilValue;

  vector<vector<VectorXd>> trajs_chains;

  auto timer = std::chrono::steady_clock::now();

  Rcpp::List outputs;
#ifdef _OPENMP
  const int burnin = control_opt["burnin"];

  int n_slope_check = (control_opt["n_slope_check"]);
  double std_lim = (control_opt["std_lim"]);
  double trend_lim = (control_opt["trend_lim"]);

  int n_chains = (control_opt["n_parallel_chain"]);
  int n_batch = (control_opt["stop_points"]);
  double start_sd = (control_opt["start_sd"]);
  double print_check_info = (control_opt["print_check_info"]);
  double max_R_hat = (control_opt["max_R_hat"]);

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

      if (n_slope_check <= curr_batch + 1)
        converge = check_conv(means, vars, curr_batch, n_slope_check, std_lim,
                              trend_lim, par_names, print_check_info,
                              batch_steps, max_R_hat, R_hat);
      all_converge =
          std::find(begin(converge), end(converge), false) == end(converge);

      // 2. if some parameter converge, stop compute gradient, or slow down the
      // gradient. if (auto_stop)
      //     for (int i=0; i < n_chains; i++) {
      //         ngmes[i]->check_converge(converge);
      //     }
      // Pflug diagnostic check for parallel chains
      if (pflug_conv_check) {
        bool pflug_converged = true;
        for (int i = 0; i < n_chains; i++) {
          // Check if max_pflug_sum has increased
          if (opt_vec[i].get_max_pflug_sum() > max_pflug_sum_before[i]) {
            pflug_converged = false;
            break;
          }
        }
        if (pflug_converged) {
          all_converge = true;
          if (verbose_enabled)
            Rcpp::Rcout << "Pflug diagnostic satisfied: max_pflug_sum did not "
                           "increase in this batch for all chains."
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
    outputs.push_back(ngmes[i]->output());
    if (store_traj) {
      opt_vec[i].record_current_state();
      trajs_chains.push_back(opt_vec[i].get_trajs());
    }
  }
  if (all_converge)
    std::cout << "Reach convergence in " << steps << " iterations."
              << std::endl;

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
std::vector<bool> check_conv(const MatrixXd &means, const MatrixXd &vars,
                             int curr_batch, int n_slope_check, double std_lim,
                             double trend_lim,
                             const std::vector<std::string> &par_names,
                             bool print_check_info, int batch_steps,
                             double max_R_hat, const VectorXd &R_hat) {
  int n_params = means.cols();
  std::vector<bool> conv(n_params, true);

  if (print_check_info) {
    std::cout << "\nstop " << curr_batch + 1 << ":\n";

    // Calculate dynamic line width
    // "Param:   " is 9 chars
    // Each param is setw(9) + 1 space = 10 chars
    int line_width = 9 + n_params * 10;

    // Print separator line
    std::cout << std::string(line_width, '-') << "\n";

    // Print parameter names
    std::cout << "Param:   ";
    for (const auto &name : par_names) {
      std::cout << std::setw(9) << std::left << name << " ";
    }
    std::cout << "\n";

    // Print separator line
    std::cout << std::string(line_width, '-') << "\n";

    // Print R-hat values with proper alignment
    std::cout << "R_hat:   ";
    for (int i = 0; i < n_params; i++) {
      std::cout << std::setw(9) << std::fixed << std::setprecision(3)
                << std::left << R_hat(i) << " ";
    }
    std::cout << "\n";

    // Print separator line
    std::cout << std::string(line_width, '-') << "\n\n";
  }

  for (int i = 0; i < n_params; i++) {
    if (R_hat(i) > max_R_hat) {
      conv[i] = false;
    }
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
