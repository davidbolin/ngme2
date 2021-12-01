#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "block.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List predict_cpp(Rcpp::List in_list) {
    // *****************   Read From Input   *****************  
    
    //observations and latents
    Rcpp::List obs_list     = Rcpp::as<Rcpp::List> (in_list["observe_in"]);
    Eigen::VectorXd Y       = obs_list["Y"];

    Rcpp::List latents_list = Rcpp::as<Rcpp::List> (in_list["latents_in"]);
    
    // config_list
    Rcpp::List config_list  = Rcpp::as<Rcpp::List> (in_list["config_in"]);

    // Flag to Specify what parameter to optimize
    const int OPT_m = config_list["opt_m"];
    const int OPT_K = config_list["opt_k"];
    const int OPT_V = config_list["opt_v"];

    BlockModel block (Y, latents_list);

    // *****************   Main Process - Optimization *****************  
    
    Optimizer opt;

    // opt.sgd(block, 0.01, 0.01, false);

    // *****************   Construct Output   *****************  
    Rcpp::List out_list;
    out_list = block.testResult();

    return out_list;
}
