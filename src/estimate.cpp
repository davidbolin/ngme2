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
                             double, double, std::string, bool, int, double,
                             const VectorXd &);

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

  Rcpp::List output = R_NilValue;

  vector<vector<VectorXd>> trajs_chains;

  auto timer = std::chrono::steady_clock::now();

  Rcpp::List outputs;
#ifdef _OPENMP
  const int burnin = control_opt["burnin"];
  const bool exchange_VW = control_opt["exchange_VW"];

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
    if (verbose_enabled && i > 0) {
      opt_vec.back().set_verbose(false);
    }
  }

  std::string par_string = ngmes[0]->get_par_string();

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

  std::vector<bool> converge(n_params, false);
  bool all_converge = false;
  int steps = 0;
  int batch_steps = (iterations / n_batch);

  int curr_batch = 0;

  while (steps < iterations && !all_converge) {
    MatrixXd mat(n_chains, n_params);
    batch_sum.setZero();
    batch_sq_sum.setZero();

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
      // exchange VW
      // if (exchange_VW) {
      //     vector<vector<VectorXd>> tmp = ngmes[0]->get_VW();
      //     for (int i = 0; i < n_chains - 1; i++) {
      //         vector<vector<VectorXd>> VW = ngmes[i+1]->get_VW();
      //         ngmes[i]->set_prev_VW(VW);
      //     }
      //     ngmes[n_chains - 1]->set_prev_VW(tmp);
      // }

      // 2. convergence            // Compute R_hat
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

      if (n_slope_check <= curr_batch + 1)
        converge = check_conv(means, vars, curr_batch, n_slope_check, std_lim,
                              trend_lim, par_string, print_check_info,
                              batch_steps, max_R_hat, R_hat);
      all_converge =
          std::find(begin(converge), end(converge), false) == end(converge);

      // 3. if some parameter converge, stop compute gradient, or slow down the
      // gradient. if (auto_stop)
      //     for (int i=0; i < n_chains; i++) {
      //         ngmes[i]->check_converge(converge);
      //     }
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
      trajs_chains.push_back(opt_vec[i].get_trajs());
    }
  }
  if (all_converge)
    std::cout << "Reach convergence in " << steps << " iterations."
              << std::endl;

#else // No parallel chain
  Ngme ngme(R_ngme, seed, sampling_strategy);
  Ngme_optimizer opt(control_opt, std::make_shared<Ngme>(ngme));
  opt.sgd(0.1, iterations, max_relative_step, max_absolute_step,
          compute_precond_each_iter);
  // estimation done, posterior sampling
  // ngme.sampling(10, true);
  outputs.push_back(ngme.output());
  if (store_traj) {
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
  return outputs;
}

// A is the horizontal concat version of latent A
// Rcpp::List sampling_cpp_given_A(const Rcpp::List& ngme_replicate, int n, bool
// posterior, unsigned long seed, const Eigen::SparseMatrix<double>& A) {
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
                             double trend_lim, std::string par_string,
                             bool print_check_info, int batch_steps,
                             double max_R_hat, const VectorXd &R_hat) {
  int n_params = means.cols();
  std::vector<bool> conv(n_params, true);

  // std::string std_line = "  < std_lim:  ";
  // std::string trend_line = " < trend_lim: ";
  std::string r_hat_line = " R_hat:        ";

  // 0. Gelman-Rubin statistic
  for (int i = 0; i < n_params; i++) {
    if (R_hat(i) > max_R_hat) {
      conv[i] = false;
    }
    char buffer[50];
    sprintf(buffer, "   %.3f", R_hat(i));
    r_hat_line += buffer;
  }
  // 1. check coef. of var. of every parameter < std_lim
  // for (int i = 0; i < n_params; i++) {
  //   // std::cout << "i = " << i << std::endl;
  //   // std::cout << "std_lim ratio = " << sqrt(vars(curr_batch, i)) /
  //   // (abs(means(curr_batch, i)) + pow(10,-5)) / std_lim << std::endl;
  //   if (sqrt(vars(curr_batch, i)) / (abs(means(curr_batch, i)) + pow(10, -5))
  //   >
  //       std_lim) {
  //     conv[i] = false;
  //     std_line += "   false"; // of length 8
  //   } else {
  //     std_line += "    true";
  //   }
  // }

  // 2. check the slope of every para < threshold
  // MatrixXd B(n_slope_check, 2); // B is the design matrix
  // B.col(0) = VectorXd::Ones(n_slope_check);
  // for (int i = 0; i < n_slope_check; i++)
  //   B(i, 1) = i;
  // // B(i, 1) = i * batch_steps;

  // for (int i = 0; i < n_params; i++) {
  //   // VectorXd mean = means.col(i)(Eigen::seq(curr_batch - n_slope_check +
  //   1,
  //   // curr_batch)); // Eigen 3.4 Eigen::seq
  //   VectorXd mean = means.block(curr_batch - n_slope_check + 1, i,
  //                               n_slope_check, 1); // Eigen block API
  //   VectorXd Sigma_inv =
  //       vars.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1)
  //           .cwiseInverse();

  //   MatrixXd Q = B.transpose() * Sigma_inv.asDiagonal() * B;
  //   Vector2d beta =
  //       Q.llt().solve(B.transpose() * Sigma_inv.asDiagonal() * mean);
  //   // MatrixXd Q = B.transpose() * B;
  //   // Vector2d beta = Q.llt().solve(B.transpose() * mean);

  //   // check the criterion
  //   // std::cout << "i = " << i << std::endl;
  //   // std::cout << "beta(1)/trend = " << abs(beta(1))/trend_lim <<
  //   std::endl;
  //   // std::cout << "Q = \n" << Q << std::endl;
  //   // std::cout << "-------" << std::endl;

  //   // Q(1,1) is way too big
  //   // if (abs(beta(1)) - 2 * sqrt(Q(1, 1)) > trend_lim * abs(beta(0))) {

  //   if (abs(beta(1)) > trend_lim) {
  //     conv[i] = false;
  //     trend_line += "   false";
  //   } else {
  //     trend_line += "    true";
  //   }
  //   // std::cout << "mean here = " << mean << std::endl;
  //   // std::cout << "Q here = " << Q << std::endl;
  //   // std::cout << "beta here = " << Q << std::endl;
  // }

  if (print_check_info)
    std::cout << "stop " << curr_batch + 1 << ": \n"
              << par_string
              << "\n"
              // << std_line << "\n"
              // << trend_line << "\n"
              << r_hat_line << "\n\n";
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
//         slope.satisfied <- abs(beta[2])-2*sqrt(Q[2,2])<trend.lim*abs(beta[1])
//         #no significant trend output[i] = std.satisfied&slope.satisfied
//       }
//     }
//     return(output)
