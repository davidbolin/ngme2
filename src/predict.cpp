#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"
#include "include/timer.h"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::SparseMatrix;

using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List predict_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list         = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const Eigen::MatrixXd X = Rcpp::as<MatrixXd>   (gen_list["X"]);
// std::cout << X << std::endl;
        const int n_regs        = Rcpp::as<int>   (gen_list["n_regs"]);
        const string family     = Rcpp::as<string>   (gen_list["family"]);
        const bool opt_fix_effect = Rcpp::as<bool>   (gen_list["opt_fix_effect"]);

    // init list
    Rcpp::List inits         = Rcpp::as<Rcpp::List> (gen_list["init"]);
        VectorXd beta = Rcpp::as<VectorXd>   (inits["beta"]);
        double sigma_eps = Rcpp::as<double>   (inits["sigma_eps"]);
    
    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);
        // Flag to Specify what parameter to optimize
        const int iterations = config_list["iterations"];
        const int n_gibbs = config_list["gibbs_sample"];
        const double stepsize = config_list["stepsize"];
        
        // const Eigen::VectorXd trueV = Rcpp::as<VectorXd>   (config_list["trueV"]);
        // const Eigen::VectorXd trueW = Rcpp::as<VectorXd>   (config_list["trueW"]);
    BlockModel block (X, Y, family, n_regs, latents_list, n_gibbs, beta, sigma_eps, opt_fix_effect);

    // *****************   Main Process - Optimization *****************  
    Optimizer opt;


    // *****************   Construct Output   ***************** 
auto timer = std::chrono::steady_clock::now();
    Rcpp::List out_list = opt.sgd(block, stepsize, 0.1, false, iterations);
std::cout << "total time is (ms): " << since(timer).count() << std::endl;   

// Rcpp::List out_list;
    return out_list;
}

