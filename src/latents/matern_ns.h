#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;
// get_K_params, grad_K_params, set_K_params, output
class matern_ns : public Latent {
public:
    matern_ns(Rcpp::List latent_in) 
    : Latent(latent_in)
    {
if (debug) std::cout << "constructor of matern ns" << std::endl;
        Rcpp::List operator_in = Rcpp::as<Rcpp::List> (latent_in["operator_in"]); // containing C and G
        ope = new nonstationaryGC(operator_in);
        
        // Init K and Q
        SparseMatrix<double> K = getK();
        SparseMatrix<double> Q = K.transpose() * K;
        
        solver_K.init(n_mesh, 0,0,0);
        solver_K.analyze(K);
        // compute_trace();

        // Init Q
        solver_Q.init(n_mesh, 0,0,0);
        solver_Q.analyze(Q);

        // fix sigma 
        sigma = 1;
        opt_flag[2] = false; 

if (debug) std::cout << "finish constructor of matern ns" << std::endl;
    }
    
    // inherit get_K_parameter, grad_K_parameter, set_K_parameter

    // generating output
    Rcpp::List get_estimates() const {
        return Rcpp::List::create(
            Rcpp::Named("theta in Matern operator") = ope->get_parameter(),
            Rcpp::Named("mu")    = mu,
            Rcpp::Named("sigma") = sigma,
            Rcpp::Named("var")   = var->get_var()
        );
    }
};

