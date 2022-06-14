#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"
#include "include/timer.h"

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
    Rcpp::List general_in    = Rcpp::as<Rcpp::List> (in_list["general_in"]);
    Rcpp::List init_values   = Rcpp::as<Rcpp::List> (in_list["init_values"]);
    Rcpp::List latents_list  = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    Rcpp::List control_list  = Rcpp::as<Rcpp::List> (in_list["control_in"]);
        const int iterations = control_list["iterations"];
    Rcpp::List debug_list  = Rcpp::as<Rcpp::List> (in_list["debug"]);
    
    BlockModel block (
        general_in, 
        latents_list, 
        control_list,
        init_values, 
        debug_list);

    // *****************   Main Process - Optimization *****************  
    Optimizer opt;
    
auto timer = std::chrono::steady_clock::now();
    // Rcpp::List trajectory = opt.sgd(block, stepsize, 0.1, false, iterations);
    Rcpp::List trajectory = opt.sgd(block, 0.1, iterations);
std::cout << "Total time is (ms): " << since(timer).count() << std::endl;   

    // *****************   Construct Output   ***************** 
    return Rcpp::List::create(
        Rcpp::Named("trajectory") = trajectory,
        // final estimate from block model
        Rcpp::Named("estimates") = block.get_estimates()
    );
}

