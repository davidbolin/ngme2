#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block_reps.h"
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
Rcpp::List estimate_cpp(const Rcpp::List& list_ngmes, const Rcpp::List& control_opt) {
    unsigned long seed = Rcpp::as<unsigned long> (control_opt["seed"]);
    std::mt19937 rng (seed);

    const int iterations = control_opt["iterations"];
    const double max_relative_step = control_opt["max_relative_step"];
    const double max_absolute_step = control_opt["max_absolute_step"];

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
    double print_check_info = (control_opt["print_check_info"]);
    omp_set_num_threads(n_chains);

    // init each model
    std::vector<std::unique_ptr<Block_reps>> block_reps;
    int i = 0;
    for (i=0; i < n_chains; i++) {
        block_reps.push_back(std::make_unique<Block_reps>(list_ngmes, rng()));
    }
    std::string par_string = block_reps[0]->get_par_string();

    // burn in period
    // #pragma omp parallel for schedule(static)
    // for (i=0; i < n_chains; i++)
        // (block_reps[i])->burn_in(burnin+3);
    int n_params = block_reps[0]->get_n_params();
    MatrixXd means (n_batch, n_params);
    MatrixXd vars (n_batch, n_params);

    std::vector<bool> converge (n_params, false);
    bool all_converge = false;
    int steps = 0;
    int batch_steps = (iterations / n_batch);

    int curr_batch = 0;

    // Not thread-safe using Rcpp::List to init optimizer
    // so better init first
    vector<Optimizer> opt_vec;
    for (int i = 0; i < n_chains; i++)
        opt_vec.push_back(Optimizer(control_opt));

    while (steps < iterations && !all_converge) {
        MatrixXd mat (n_chains, n_params);
        #pragma omp parallel for schedule(static)
        for (i=0; i < n_chains; i++) {
            // Optimizer opt (control_opt);
            // VectorXd param = opt.sgd(*(blocks[i]), 0.1, batch_steps, max_relative_step, max_absolute_step);
            VectorXd param = opt_vec[i].sgd(*(block_reps[i]), 0.1, batch_steps, max_relative_step, max_absolute_step);

            #pragma omp critical
            mat.row(i) = param;
        }
        steps += batch_steps;

        // compute mean and variance
        means.row(curr_batch) = mat.colwise().mean();
        for (int k=0; k < n_params; k++)
            vars(curr_batch, k) = (mat.col(k).array() - means(curr_batch, k)).square().sum() / (n_chains - 1);

        if (n_chains > 1) {
            // exchange VW
            if (exchange_VW) {
                vector<vector<VectorXd>> tmp = block_reps[0]->get_VW();
                for (int i = 0; i < n_chains - 1; i++) {
                    vector<vector<VectorXd>> VW = block_reps[i+1]->get_VW();
                    block_reps[i]->set_prev_VW(VW);
                }
                block_reps[n_chains - 1]->set_prev_VW(tmp);
            }

            // 2. convergence check
            if (n_slope_check <= curr_batch + 1)
                converge = check_conv(means, vars, curr_batch, n_slope_check, std_lim, trend_lim, par_string, print_check_info);
            all_converge = std::find(begin(converge), end(converge), false) == end(converge);

            // 3. if some parameter converge, stop compute gradient, or slow down the gradient.
            // if (auto_stop)
            //     for (int i=0; i < n_chains; i++) {
            //         block_reps[i]->check_converge(converge);
            //     }
        }

        curr_batch++;
    }

// After estimation
// ****** posterior sampling (sampling each chain..)
    // #pragma omp parallel for schedule(static)
    // for (i=0; i < n_chains; i++) {
    //     block_reps[i]->sampling(100, true);
    // }

    // generate outputs
    for (i=0; i < n_chains; i++) {
        outputs.push_back(block_reps[i]->output());
        trajs_chains.push_back(opt_vec[i].get_trajs());
    }
    if (all_converge)
        std::cout << "Reach convergence in " << steps << " iterations." << std::endl;
    else
        std::cout << "Estimation ends." << std::endl;

#else // No parallel chain
    Block_reps block_reps (list_ngmes, rng());
    Optimizer opt (control_opt);
    opt.sgd(block_reps, 0.1, iterations, max_relative_step, max_absolute_step);
    // estimation done, posterior sampling
    // block_reps.sampling(10, true);
    outputs.push_back(block_reps.output());
    trajs_chains.push_back(opt.get_trajs());
#endif

std::cout << "Total time of the estimation is (s): " << since(timer).count() / 1000 << std::endl;

    outputs.attr("opt_traj") = trajs_chains;
    return outputs;
}

// [[Rcpp::export]]
Rcpp::List sampling_cpp(const Rcpp::List& ngme_block, int n, bool posterior, unsigned long seed) {
    std::mt19937 rng (seed);
    BlockModel block (ngme_block, rng());

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