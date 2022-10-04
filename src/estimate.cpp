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
#include <random>

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List estimate_cpp(Rcpp::List ngme_block) {
    unsigned long seed = Rcpp::as<unsigned long> (ngme_block["seed"]);
    std::mt19937 rng (seed);

    Rcpp::List control_in = ngme_block["control"];
    const int iterations = (control_in["iterations"]);
    int n_chains = (control_in["n_parallel_chain"]);

    Rcpp::List trajectory = R_NilValue;
    Rcpp::List output = R_NilValue;

    Rcpp::List out;
auto timer = std::chrono::steady_clock::now();

#ifdef _OPENMP
    // Nb break points
    int Nb = 10;
    omp_set_num_threads(n_chains);
    #pragma omp parallel
    {
        // unsigned long chain_seed = rng();
        int i = omp_get_thread_num();
        // unsigned long chain_seed = rng();
        BlockModel block (ngme_block, i);

        Optimizer opt;
        trajectory = opt.sgd(block, 0.1, iterations);
            // ngme.attr("trajectory") = trajectory;
        #pragma omp critical
            out[(String) (i+1)] = block.output();
    }
#else
    BlockModel block (ngme_block, rng());
    Optimizer opt;
    trajectory = opt.sgd(block, 0.1, iterations);
    Rcpp::List ngme = block.output();
    // ngme.attr("trajectory") = trajectory;
    out[1] = ngme;
#endif

std::cout << "Total time is (ms): " << since(timer).count() << std::endl;

    // return Rcpp::List::create(
    //     Rcpp::Named("opt_trajectory") = trajectory,
    //     Rcpp::Named("estimation") = output
    // );
    return out;
}

// [[Rcpp::export]]
Rcpp::List sampling_cpp(Rcpp::List& ngme_block, int iterations, bool posterior) {
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