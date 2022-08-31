#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"
#include "timer.h"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List estimate_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    Rcpp::List general_list  = Rcpp::as<Rcpp::List> (in_list["general_in"]);
    Rcpp::List latents_list  = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    Rcpp::List noise_list    = Rcpp::as<Rcpp::List> (in_list["noise_in"]);
    Rcpp::List control_list  = Rcpp::as<Rcpp::List> (in_list["control_in"]);
        const int iterations = control_list["iterations"];
    Rcpp::List debug_list    = Rcpp::as<Rcpp::List> (in_list["debug"]);
    
    BlockModel block (
        general_list, 
        latents_list, 
        noise_list, 
        control_list,
        debug_list
    );

    Rcpp::List trajectory = R_NilValue;
    bool estimation = Rcpp::as<bool>    (control_list["estimation"]);
    if (!estimation) { 
        // no optimization, only gibbs sampling
        int steps  = Rcpp::as<int> (control_list["burnin"]) + Rcpp::as<int> (control_list["iterations"]);
        block.burn_in(steps);
    } else {    
        // doing optimization
auto timer = std::chrono::steady_clock::now();
        Optimizer opt;
        Rcpp::List trajectory = opt.sgd(block, 0.1, iterations);
std::cout << "Total time is (ms): " << since(timer).count() << std::endl;
    }

    return Rcpp::List::create(
            Rcpp::Named("trajectory") = trajectory,
            Rcpp::Named("output") = block.output()
    );
}

