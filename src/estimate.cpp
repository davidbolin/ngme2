#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "ngme.h"
#include "include/timer.h"

#ifdef _OPENMP
    #include<omp.h>
#endif

#include <iostream>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <random>

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using namespace Rcpp;

std::vector<bool> check_conv(const MatrixXd&, const MatrixXd&, int, int, double, double, std::string, bool);

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List estimate_cpp(const Rcpp::List& R_ngme, const Rcpp::List& control_opt) {
    unsigned long seed = Rcpp::as<unsigned long> (control_opt["seed"]);
    std::mt19937 rng (seed);

    const int iterations = control_opt["iterations"];
    const double max_relative_step = control_opt["max_relative_step"];
    const double max_absolute_step = control_opt["max_absolute_step"];
    const int sampling_strategy = control_opt["sampling_strategy"];
    bool compute_precond_each_iter = true;

    Rcpp::List output = R_NilValue;

    vector<vector<VectorXd>> trajs_chains;

auto timer = std::chrono::steady_clock::now();

    Rcpp::List outputs;
#ifdef _OPENMP
    const int burnin = control_opt["burnin"];
    const bool exchange_VW = control_opt["exchange_VW"];
    const bool precond_by_diff_chain = control_opt["precond_by_diff_chain"];

    int n_slope_check = (control_opt["n_slope_check"]);
    double std_lim = (control_opt["std_lim"]);
    double trend_lim = (control_opt["trend_lim"]);

    int n_chains = (control_opt["n_parallel_chain"]);
    int n_batch = (control_opt["stop_points"]);
    double print_check_info = (control_opt["print_check_info"]);
    omp_set_num_threads(n_chains);

    if (n_chains > 1 && precond_by_diff_chain) {
        compute_precond_each_iter = false;
    }

    // init model and optimizer
    vector<std::shared_ptr<Ngme>> ngmes;
    vector<Ngme_optimizer> opt_vec;
    int i = 0;
    for (i=0; i < n_chains; i++) {
        // Not thread-safe using Rcpp::List to init optimizer
        ngmes.push_back(std::make_shared<Ngme>(R_ngme, rng(), sampling_strategy));
        opt_vec.push_back(Ngme_optimizer(control_opt, ngmes[i]));
    }

    std::string par_string = ngmes[0]->get_par_string();

    // burn in period
    // #pragma omp parallel for schedule(static)
    // for (i=0; i < n_chains; i++)
        // (ngmes[i])->burn_in(burnin+3);
    int n_params = ngmes[0]->get_n_params();
    MatrixXd means (n_batch, n_params);
    MatrixXd vars (n_batch, n_params);

    std::vector<bool> converge (n_params, false);
    bool all_converge = false;
    int steps = 0;
    int batch_steps = (iterations / n_batch);

    int curr_batch = 0;

    while (steps < iterations && !all_converge) {
        MatrixXd mat (n_chains, n_params);
        #pragma omp parallel for schedule(static)
        for (i=0; i < n_chains; i++) {
            VectorXd param = opt_vec[i].sgd(
                0.1,
                batch_steps,
                max_relative_step,
                max_absolute_step,
                compute_precond_each_iter
            );

// compute preconditoner
// update A1.precond = 1/3 (A2.precond + A3.precond + A4.precond)
            #pragma omp critical
            mat.row(i) = param;
        }
        steps += batch_steps;

        // compute mean and variance
        means.row(curr_batch) = mat.colwise().mean();
        for (int k=0; k < n_params; k++)
            vars(curr_batch, k) = (mat.col(k).array() - means(curr_batch, k)).square().sum() / (n_chains - 1);

        if (n_chains > 1) {
            if (precond_by_diff_chain) {
                // set the preconditioner by other chains
                MatrixXd precond_sum = MatrixXd::Zero(n_params, n_params);
                for (int i=0; i < n_chains; i++)
                    precond_sum += opt_vec[i].get_preconditioner();
                for (int i=0; i < n_chains; i++)
                    opt_vec[i].set_preconditioner((precond_sum - opt_vec[i].get_preconditioner()) / (n_chains - 1));
            }

            // exchange VW
            if (exchange_VW) {
                vector<vector<VectorXd>> tmp = ngmes[0]->get_VW();
                for (int i = 0; i < n_chains - 1; i++) {
                    vector<vector<VectorXd>> VW = ngmes[i+1]->get_VW();
                    ngmes[i]->set_prev_VW(VW);
                }
                ngmes[n_chains - 1]->set_prev_VW(tmp);
            }

            // 2. convergence check
            if (n_slope_check <= curr_batch + 1)
                converge = check_conv(means, vars, curr_batch, n_slope_check, std_lim, trend_lim, par_string, print_check_info);
            all_converge = std::find(begin(converge), end(converge), false) == end(converge);

            // 3. if some parameter converge, stop compute gradient, or slow down the gradient.
            // if (auto_stop)
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

    // generate outputs
    for (i=0; i < n_chains; i++) {
        outputs.push_back(ngmes[i]->output());
        trajs_chains.push_back(opt_vec[i].get_trajs());
    }
    if (all_converge)
        std::cout << "Reach convergence in " << steps << " iterations." << std::endl;

#else // No parallel chain
    Ngme ngme (R_ngme, rng(), sampling_strategy);
    Ngme_optimizer opt (control_opt, std::make_shared<Ngme>(ngme));
    opt.sgd(
        0.1,
        iterations,
        max_relative_step, max_absolute_step,
        compute_precond_each_iter
    );
    // estimation done, posterior sampling
    // ngme.sampling(10, true);
    outputs.push_back(ngme.output());
    trajs_chains.push_back(opt.get_trajs());
#endif

std::cout << "Estimation ends." << std::endl;
std::cout << "Total time of the estimation is (s): " << since(timer).count() / 1000 << std::endl;

    outputs.attr("opt_traj") = trajs_chains;
    return outputs;
}

// [[Rcpp::export]]
Rcpp::List sampling_cpp(const Rcpp::List& ngme_replicate, int n, bool posterior, unsigned long seed) {
    std::mt19937 rng (seed);
    BlockModel block (ngme_replicate, rng());

    return block.sampling(n, posterior);
}

/*
    For checking convergence of parallel chains
    data is n_iters () * n_params (how many params in total)
*/
std::vector<bool> check_conv(
    const MatrixXd& means,
    const MatrixXd& vars,
    int curr_batch,
    int n_slope_check,
    double std_lim,
    double trend_lim,
    std::string par_string,
    bool print_check_info
) {
    int n_params = means.cols();
    std::vector<bool> conv (n_params, true);

    std::string std_line    = "  < std_lim:  ";
    std::string trend_line  = " < trend_lim: ";

    // 1. check coef. of var. of every parameter < std_lim
    for (int i=0; i < n_params; i++)
        if (sqrt(vars(curr_batch, i)) / (abs(means(curr_batch, i)) + pow(10,-5)) > std_lim) {
            conv[i] = false;
            std_line += "   false"; // of length 8
        } else {
            std_line += "    true";
        }

    // 2. check the slope of every para < threshold
    MatrixXd B (n_slope_check, 2);
        B.col(0) = VectorXd::Ones(n_slope_check);
        for (int i=0; i < n_slope_check; i++)
            B(i, 1) = i;
    for (int i = 0; i < n_params; i++) {
        // VectorXd mean = means.col(i)(Eigen::seq(curr_batch - n_slope_check + 1, curr_batch)); // Eigen 3.4 Eigen::seq
        VectorXd mean      = means.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1);  // Eigen block API
        VectorXd Sigma_inv = vars.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1).cwiseInverse();
        MatrixXd Q = B.transpose() * Sigma_inv.asDiagonal() * B;
        Vector2d beta = Q.llt().solve(B.transpose() * Sigma_inv.asDiagonal() * mean);
        if (abs(beta(1)) - 2 * sqrt(Q(1, 1)) > trend_lim * abs(beta(0))) {
            conv[i] = false;
            trend_line += "   false";
        } else {
            trend_line += "    true";
        }
// std::cout << "mean here = " << mean << std::endl;
// std::cout << "Q here = " << Q << std::endl;
// std::cout << "beta here = " << Q << std::endl;
    }

    if (print_check_info) std::cout << "stop " << curr_batch+1 << ": \n"
        << par_string << "\n"
        << std_line << "\n"
        << trend_line << "\n\n";
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
//         slope.satisfied <- abs(beta[2])-2*sqrt(Q[2,2])<trend.lim*abs(beta[1]) #no significant trend
//         output[i] = std.satisfied&slope.satisfied
//       }
//     }
//     return(output)