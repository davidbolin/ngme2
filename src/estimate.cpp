#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"
#include "timer.h"

#ifdef _OPENMP
#include<omp.h>
#endif

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List estimate_cpp(Rcpp::List& ngme_block) {

    Rcpp::List control_in = ngme_block["control"];
    const int iterations = (control_in["iterations"]);
    int n_chains = (control_in["n_parallel_chain"]);

    Rcpp::List trajectory = R_NilValue;
    Rcpp::List output = R_NilValue;

    Rcpp::List result;

auto timer = std::chrono::steady_clock::now();
#ifdef _OPENMP
    omp_set_num_threads(n_chains);
    #pragma omp parallel for
    {
        for (int i=0; i < n_chains; i++) {
            BlockModel block (ngme_block);
            Optimizer opt;
            trajectory = opt.sgd(block, 0.1, iterations);
            result[(String) i] = trajectory;
        }
    }
#else
    BlockModel block (ngme_block);
    Optimizer opt;
    trajectory = opt.sgd(block, 0.1, iterations);
    output = block.output();
#endif

std::cout << "Total time is (ms): " << since(timer).count() << std::endl;


    // return Rcpp::List::create(
    //     Rcpp::Named("opt_trajectory") = trajectory,
    //     Rcpp::Named("estimation") = output
    // );
    return result;
}

// [[Rcpp::export]]
Rcpp::List sampling_cpp(Rcpp::List& ngme_block, int iterations, bool posterior) {
    BlockModel block (ngme_block);

    return block.sampling(iterations, posterior);
}

// void update_estimation(Rcpp::List& ngme_block, const BlockModel& block);
// void update_estimation(Rcpp::List& ngme_block, const BlockModel& block) {
//     // design the input of R
//     // ngme_block["V"] <- ...
// }