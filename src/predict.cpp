#include <Rcpp.h>
#include <RcppEigen.h>
#include "optimizer.h"
#include "model.h"
#include "block.h"

using namespace Rcpp;

// [[Rcpp::export]]
List predict_cpp(Rcpp::List in_list) {
    // read from in_list
    Rcpp::List observe_in   = Rcpp::as<Rcpp::List> (in_list["observe_in"]);
    Rcpp::List operator_in  = Rcpp::as<Rcpp::List> (in_list["operator_in"]);
    Rcpp::List var_in       = Rcpp::as<Rcpp::List> (in_list["var_in"]);

    // Block block;
    // Operator operator(operator_in)
    //

    
    // *****************   Main Process   *****************  
    SomeFun f;
    Optimizer opt;
    VectorXd out = opt.sgd(f, 0.05, 0.01, false);
    







    // construct output
    Rcpp::List out_list;

    out_list["operator_in"] = operator_in;
    return out_list;
}
