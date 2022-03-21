#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"
#include "include/timer.h"

#include <Eigen/Sparse>

using Eigen::SparseMatrix;

using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List predict_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list     = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y       = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const int n_regs        = Rcpp::as<int>   (gen_list["n_regs"]);

    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int interations = config_list["iterations"];
        const int n_gibbs = config_list["gibbs_sample"];
        const double stepsize = config_list["stepsize"];
        
        const Eigen::VectorXd trueV = Rcpp::as<VectorXd>   (config_list["trueV"]);
        const Eigen::VectorXd trueW = Rcpp::as<VectorXd>   (config_list["trueW"]);
    

    BlockModel block (Y, n_regs, latents_list, n_gibbs);


    // *****************   Main Process - Optimization *****************  
    
    Optimizer opt;

    // opt.sgd(block, 0.01, 0.01, false);

    // *****************   Construct Output   ***************** 


    Rcpp::List out_list;


    // *****************   testing  ***************** 

    
    // out_list = block.testResult();
    return out_list;
}


// [[Rcpp::export]]
Rcpp::List test_output(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list     = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y       = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const int n_regs        = Rcpp::as<int>   (gen_list["n_regs"]);

    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int iterations = config_list["iterations"];
        const int n_gibbs = config_list["gibbs_sample"];
        const double stepsize = config_list["stepsize"];
        
        const Eigen::VectorXd trueV = Rcpp::as<VectorXd>   (config_list["trueV"]);
        const Eigen::VectorXd trueW = Rcpp::as<VectorXd>   (config_list["trueW"]);
    BlockModel block (Y, n_regs, latents_list, n_gibbs);

    // *****************   Main Process - Optimization *****************  
    
    Optimizer opt;


    // *****************   Construct Output   ***************** 

    Rcpp::List out_list = opt.sgd(block, stepsize, 0.1, false, iterations);

    return out_list;
}

