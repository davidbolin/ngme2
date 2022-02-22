#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"

#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::Map;

using namespace Rcpp;

void testGradConvergence(BlockModel&, int);

// [[Rcpp::export]]
Rcpp::List predict_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list     = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y       = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const int n_paras       = Rcpp::as<int>   (gen_list["n_paras"]);
        const int n_regs        = Rcpp::as<int>   (gen_list["n_regs"]);

    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int OPT_m = config_list["opt_m"];
        const int OPT_K = config_list["opt_k"];
        const int OPT_V = config_list["opt_v"];
        const int interations = config_list["iterations"];
        
        const Eigen::VectorXd trueV = Rcpp::as<VectorXd>   (config_list["trueV"]);
        const Eigen::VectorXd trueW = Rcpp::as<VectorXd>   (config_list["trueW"]);
    
    BlockModel block (Y, n_paras, n_regs, latents_list);
        block.setW(trueW);

    // *****************   Main Process - Optimization *****************  
    
    Optimizer opt;

    // opt.sgd(block, 0.01, 0.01, false);

    // *****************   Construct Output   ***************** 


    Rcpp::List out_list;


    // *****************   testing  ***************** 

    // 1. test convergence
    // block.setVW(); // setting the trueV and trueW
    testGradConvergence(block, interations);
    
    // out_list = block.testResult();
    return out_list;
}

void testGradConvergence(BlockModel& block, int iterations) {
    Optimizer opt;
    opt.sgd(block, 0.05, 0.1, false, iterations);
}




// [[Rcpp::export]]
Rcpp::List test_init(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list     = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y       = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const int n_paras       = Rcpp::as<int>   (gen_list["n_paras"]);
        const int n_regs        = Rcpp::as<int>   (gen_list["n_regs"]);

    
    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int OPT_m = config_list["opt_m"];
        const int OPT_K = config_list["opt_k"];
        const int OPT_V = config_list["opt_v"];

    BlockModel block (Y, n_paras, n_regs, latents_list);
    block.sampleW_VY();
// std::cout << "theta=" << block.get_parameter() <<std::endl;
// std::cout << "grad=" << block.grad() <<std::endl;
    // block.setVW(); 

    Rcpp::List out_list;
    out_list = block.testResult();

    return out_list;
}