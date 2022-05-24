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
Rcpp::List estimate_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List gen_list         = Rcpp::as<Rcpp::List> (in_list["general_in"]);
        const Eigen::VectorXd Y = Rcpp::as<VectorXd>   (gen_list["Y"]);
        const Eigen::MatrixXd X = Rcpp::as<MatrixXd>   (gen_list["X"]);
        const int n_regs        = Rcpp::as<int>    (gen_list["n_regs"]);
        const string family     = Rcpp::as<string> (gen_list["family"]);

    // init list
    Rcpp::List inits         = Rcpp::as<Rcpp::List> (gen_list["init"]);
        VectorXd beta = Rcpp::as<VectorXd>   (inits["beta"]);
        double sigma_eps = Rcpp::as<double>   (inits["sigma_eps"]);
    
    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // control_list
    Rcpp::List control_list  = Rcpp::as<Rcpp::List> (in_list["control_in"]);
        // Flag to Specify what parameter to optimize
        const int burnin = control_list["burnin"];
        const int iterations = control_list["iterations"];
        const int n_gibbs = control_list["gibbs_sample"];
        const double stepsize = control_list["stepsize"];
        
        const bool opt_fix_effect = Rcpp::as<bool>   (control_list["opt_fix_effect"]);

    // debug list
    Rcpp::List debug_list  = Rcpp::as<Rcpp::List> (in_list["debug"]);
        const bool debug = Rcpp::as<bool> (debug_list["debug"]);
        const bool fixW  = Rcpp::as<bool> (debug_list["fixW"]);
        const bool fixSV = Rcpp::as<bool> (debug_list["fixSV"]);
        const bool fixSigEps = Rcpp::as<bool> (debug_list["fixSigEps"]);

        Eigen::VectorXd trueW, trueSV; 
        if (fixW)      trueW = Rcpp::as<VectorXd>   (debug_list["trueW"]);
        if (fixSV)     trueSV = Rcpp::as<VectorXd>  (debug_list["trueSV"]);
        if (fixSigEps) sigma_eps = Rcpp::as<double> (debug_list["sigEps"]);
    
    BlockModel block (X, Y, family, n_regs, latents_list, 
        n_gibbs, beta, sigma_eps, opt_fix_effect, burnin, 
        // debug
        debug, fixW, fixSV, fixSigEps, trueW, trueSV);

    // *****************   Main Process - Optimization *****************  
    Optimizer opt;


    // *****************   Construct Output   ***************** 
auto timer = std::chrono::steady_clock::now();
    Rcpp::List trajectory = opt.sgd(block, stepsize, 0.1, false, iterations);
std::cout << "total time is (ms): " << since(timer).count() << std::endl;   

    // final estimate from block model


// Rcpp::List out_list;
    return Rcpp::List::create(
        Rcpp::Named("trajectory") = trajectory,
        Rcpp::Named("estimates") = block.get_estimates()
    );
}

