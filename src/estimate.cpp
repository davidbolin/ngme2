#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"
#include "timer.h"

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

using namespace Rcpp;

bool check_conv(const MatrixXd&, const MatrixXd&, int, int, double, double);

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List estimate_cpp(const Rcpp::List& ngme_block) {
    unsigned long seed = Rcpp::as<unsigned long> (ngme_block["seed"]);
    std::mt19937 rng (seed);

    Rcpp::List control_in = ngme_block["control"];
    const int iterations = control_in["iterations"];
    const int burnin = control_in["burnin"];

    Rcpp::List trajectory = R_NilValue;
    Rcpp::List output = R_NilValue;

auto timer = std::chrono::steady_clock::now();

    Rcpp::List outputs;
#ifdef _OPENMP
    // Rcout << "run parallel n_chains chains"
    int n_slope_check = (control_in["n_slope_check"]);
    double std_lim = (control_in["std_lim"]);
    double trend_lim = (control_in["trend_lim"]);

    int n_chains = (control_in["n_parallel_chain"]);
    int n_batch = (control_in["stop_points"]);
    omp_set_num_threads(n_chains);

    // init each model
    std::vector<std::unique_ptr<BlockModel>> blocks;
    int i = 0;
    for (i=0; i < n_chains; i++) {
        blocks.push_back(std::make_unique<BlockModel>(ngme_block, rng()));
    }

    // burn in period
    #pragma omp parallel for schedule(static)
    for (i=0; i < n_chains; i++)
        (blocks[i])->burn_in(burnin+3);

    int n_params = blocks[0]->get_n_params();
    MatrixXd means (n_batch, n_params);
    MatrixXd vars (n_batch, n_params);

    bool converge = false;
    int steps = 0;
    int batch_steps = (iterations > n_batch) ? (iterations / n_batch) : 1;

    int curr_batch = 0;
    while (steps < iterations && !converge) {
        MatrixXd mat (n_chains, n_params);
        #pragma omp parallel for schedule(static)
        for (i=0; i < n_chains; i++) {
            Optimizer opt;
            VectorXd param = opt.sgd(*(blocks[i]), 0.1, batch_steps);

            #pragma omp critical
            mat.row(i) = param;
        }
        steps += batch_steps;

        // compute mean and variance
        means.row(curr_batch) = mat.colwise().mean();
        for (int k=0; k < n_params; k++)
            vars(curr_batch, k) = (mat.col(k).array() - means(curr_batch, k)).square().sum() / (n_chains - 1);

        if (n_chains > 1) {
            // 1. set hessian by setting prevV and prevW (round robin)
            vector<VectorXd> tmp = blocks[0]->get_VW();
            for (int i = 0; i < n_chains - 1; i++) {
                vector<VectorXd> VW = blocks[i+1]->get_VW();
                blocks[i]->set_prev_VW(VW);
            }
            blocks[n_chains - 1]->set_prev_VW(tmp);

            // 2. convergence check
            if (n_slope_check <= curr_batch + 1)
                converge = check_conv(means, vars, curr_batch, n_slope_check, std_lim, trend_lim);
        }
        curr_batch++;
    }

    // generate outputs
    for (i=0; i < n_chains; i++) {
        outputs.push_back(blocks[i]->output());
    }
    if (converge)
        Rcpp::Rcout << "Reach convergence in " << steps << " iterations." << std::endl;
    else
        Rcpp::Rcout << "Not sure about the convergence." << std::endl;

#else
    BlockModel block (ngme_block, rng());
    Optimizer opt;
    trajectory = opt.sgd(block, 0.1, iterations);
    Rcpp::List ngme = block.output();
    outputs.push_back(block.output());
#endif

std::cout << "Total time is (ms): " << since(timer).count() << std::endl;

    return outputs;
}

// [[Rcpp::export]]
Rcpp::List sampling_cpp(const Rcpp::List& ngme_block, int iterations, bool posterior) {
    unsigned long seed = Rcpp::as<unsigned long> (ngme_block["seed"]);
    std::mt19937 rng (seed);
    BlockModel block (ngme_block, rng());

    return block.sampling(iterations, posterior);
}

/*
    For checking convergence of parallel chains
    data is n_iters () * n_params (how many params in total)
*/
bool check_conv(
    const MatrixXd& means,
    const MatrixXd& vars,
    int curr_batch,
    int n_slope_check,
    double std_lim,
    double trend_lim
) {
    bool conv = true;
    int n_params = means.cols();

    // 1. check coef. of var. of every parameter < std_lim
    for (int i=0; i < n_params && conv; i++)
        if (sqrt(vars(curr_batch, i)) / abs(means(curr_batch, i)) > std_lim)
            conv = false;

    // 2. check the slope of every para < threshold
    MatrixXd B (n_slope_check, 2);
        B.col(0) = VectorXd::Ones(n_slope_check);
        for (int i=0; i < n_slope_check; i++)
            B(i, 1) = i;
    for (int i = 0; i < n_params && conv; i++) {
        // VectorXd mean = means.col(i)(Eigen::seq(curr_batch - n_slope_check + 1, curr_batch)); // Eigen 3.4 Eigen::seq
        VectorXd mean      = means.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1);  // Eigen block API
        VectorXd Sigma_inv = vars.block(curr_batch - n_slope_check + 1, i, n_slope_check, 1).cwiseInverse();
        MatrixXd Q = B.transpose() * Sigma_inv.asDiagonal() * B;
        Vector2d beta = Q.llt().solve(B.transpose() * Sigma_inv.asDiagonal() * mean);
        if (abs(beta(1)) - 2 * sqrt(Q(1, 1)) > trend_lim * abs(beta(0)))
            conv = false;

// std::cout << "mean here = " << mean << std::endl;
// std::cout << "Q here = " << Q << std::endl;
// std::cout << "beta here = " << Q << std::endl;
    }
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