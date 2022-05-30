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
    //observations and latents
    Rcpp::List gen_list         = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const VectorXd Y        = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const MatrixXd X        = Rcpp::as<MatrixXd>   (gen_list["X"]);
        const int n_regs        = Rcpp::as<int>    (gen_list["n_regs"]);
        const string family     = Rcpp::as<string> (gen_list["family"]);

    // init list
    Rcpp::List inits         = Rcpp::as<Rcpp::List> (gen_list["init"]);
        VectorXd beta = Rcpp::as<VectorXd>   (inits["beta"]);
        double sigma_eps = Rcpp::as<double>   (inits["sigma_eps"]);
    
    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // control_list
    Rcpp::List control_list  = Rcpp::as<Rcpp::List> (in_list["control_in"]);
        const int iterations = control_list["iterations"];

    // debug list
    Rcpp::List debug_list  = Rcpp::as<Rcpp::List> (in_list["debug"]);
    
    BlockModel block (
        gen_list, 
        inits, 
        latents_list, 
        control_list,
        debug_list);

    // *****************   Main Process - Optimization *****************  
    Optimizer opt;
    
    // *****************   Construct Output   ***************** 
auto timer = std::chrono::steady_clock::now();
    
    // Rcpp::List trajectory = opt.sgd(block, stepsize, 0.1, false, iterations);
    Rcpp::List trajectory = opt.sgd(block, 0.1, iterations);

std::cout << "total time is (ms): " << since(timer).count() << std::endl;   

    // final estimate from block model

// Rcpp::List out_list;
    return Rcpp::List::create(
        Rcpp::Named("trajectory") = trajectory,
        Rcpp::Named("estimates") = block.get_estimates()
    );
}

