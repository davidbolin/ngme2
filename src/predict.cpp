#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"

#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::Map;

using namespace Rcpp;

void testGradConvergence(BlockModel& );

// [[Rcpp::export]]
Rcpp::List predict_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list     = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y       = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const int n_paras       = Rcpp::as<int>   (gen_list["n_paras"]);
        const int n_reg         = Rcpp::as<int>   (gen_list["n_reg"]);
    
    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int OPT_m = config_list["opt_m"];
        const int OPT_K = config_list["opt_k"];
        const int OPT_V = config_list["opt_v"];

    BlockModel block (Y, n_paras, n_reg, latents_list);


    // *****************   Main Process - Optimization *****************  
    
    Optimizer opt;

    // opt.sgd(block, 0.01, 0.01, false);

    // *****************   Construct Output   ***************** 


    Rcpp::List out_list;
    out_list = block.testResult();

    // *****************   testing  ***************** 

    // 1. test convergence
    // block.testGrad(); // setting the trueV and trueW
    // testGradConvergence(block);

    return out_list;
}

void testGradConvergence(BlockModel& block) {
    Optimizer opt;
    opt.sgd(block, 0.005, 0.1, false, 100);
}




// [[Rcpp::export]]
Rcpp::List test_init(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list     = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y       = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const int n_paras       = Rcpp::as<int>   (gen_list["n_paras"]);
        const int n_reg         = Rcpp::as<int>   (gen_list["n_reg"]);

    
    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int OPT_m = config_list["opt_m"];
        const int OPT_K = config_list["opt_k"];
        const int OPT_V = config_list["opt_v"];

    BlockModel block (Y, n_paras, n_reg, latents_list);


    Rcpp::List out_list;
    out_list = block.testResult();

    return out_list;
}