#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "model.h"
#include "block.h"
#include "latent.h"
// #include "var.h"

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
    
    Block block (Y, latents_list);

    // *****************   Main Process   *****************  
    
    Optimizer opt;
    // opt.sgd(block, 0.01, 0.01, false);
    

    // *****************   Construct Output   *****************  
    Rcpp::List out_list;
    out_list = block.result();

    return out_list;
}
