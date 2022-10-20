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

bool check_conv(const std::vector<VectorXd>&);

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List estimate_cpp(const Rcpp::List& ngme_block) {
    unsigned long seed = Rcpp::as<unsigned long> (ngme_block["seed"]);
    std::mt19937 rng (seed);

    Rcpp::List control_in = ngme_block["control"];
    const int iterations = (control_in["iterations"]);

    Rcpp::List trajectory = R_NilValue;
    Rcpp::List output = R_NilValue;

auto timer = std::chrono::steady_clock::now();

    Rcpp::List outputs;
#ifdef _OPENMP
    // Rcout << "run parallel n_chains chains"
    int n_chains = (control_in["n_parallel_chain"]);
    int n = (control_in["stop_points"]);
    omp_set_num_threads(n_chains);

    std::vector<std::unique_ptr<BlockModel>> blocks;
    int i = 0;
    for (i=0; i < n_chains; i++) {
        blocks.push_back(std::make_unique<BlockModel>(ngme_block, rng()));
    }

    std::vector<VectorXd> grads (n_chains);

    bool converge = false;
    int steps = 0;

    while (steps < iterations && !converge) {
        #pragma omp parallel for schedule(static)
        for (i=0; i < n_chains; i++) {
            Optimizer opt;
            VectorXd grad = opt.sgd(*(blocks[i]), 0.1, iterations / n);
            #pragma omp critical
            grads[i] = grad;

        }
        // if (omp_get_num_thread==0) {

        // // check
        // }
        steps += iterations / n;

        // set hessian by setting prevV and prevW (round robin)
        if (n_chains > 1) {
            vector<VectorXd> tmp = blocks[0]->get_VW();
            for (int i = 0; i < n_chains - 1; i++) {
                vector<VectorXd> VW = blocks[i+1]->get_VW();
                blocks[i]->set_prev_VW(VW);
            }
            blocks[n_chains - 1]->set_prev_VW(tmp);

            // check grads for convergence (t-test)
            converge = check_conv(grads);
        }
    }
    // generate outputs
    for (i=0; i < n_chains; i++) {
        outputs.push_back(blocks[i]->output());
    }

#else
    BlockModel block (ngme_block, rng());
    Optimizer opt;
    trajectory = opt.sgd(block, 0.1, iterations);
    Rcpp::List ngme = block.output();
    outputs.push_back(block.output());
#endif

std::cout << "Total time is (ms): " << since(timer).count() << std::endl;

    // return Rcpp::List::create(
    //     Rcpp::Named("opt_trajectory") = trajectory,
    //     Rcpp::Named("estimation") = output
    // );
    return outputs;
}

// [[Rcpp::export]]
Rcpp::List sampling_cpp(const Rcpp::List& ngme_block, int iterations, bool posterior) {
    unsigned long seed = Rcpp::as<unsigned long> (ngme_block["seed"]);
    std::mt19937 rng (seed);
    BlockModel block (ngme_block, rng());

    return block.sampling(iterations, posterior);
}

// void update_estimation(Rcpp::List& ngme_block, const BlockModel& block);
// void update_estimation(Rcpp::List& ngme_block, const BlockModel& block) {
//     // design the input of R
//     // ngme_block["V"] <- ...
// }

bool check_conv(const std::vector<VectorXd>& grads) {
    bool conv = true;
    int n_chains = grads.size();
    int n = grads[1].size();

    Rcpp::NumericVector v (n_chains);

    // save the mean and sd at every break point
    // check sd < eps
    for (int j=0; j < n; j++) {
        for (int i=0; i < n_chains; i++) {
            v[i] = grads[i][j];
        }
        double t = (double) Rcpp::mean(v) / (double) Rcpp::sd(v);

        // p-value: 2 * min{cdf(x) , 1 - cdf(x)}
        double cdf = R::pt(t, n_chains-1, 1, 0);
        double p_val = (cdf < 1-cdf) ? 2*cdf : 2*(1-cdf);

        if (p_val < 0.05) conv = false;
        std::cout << "p_val = " << p_val << std::endl;
    }

    // for (int i=0; i<n_chains; i++) {
    //     std::cout << "grads[i] = " << grads[i] << std::endl;
    // }
    return false;
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