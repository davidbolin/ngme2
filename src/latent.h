#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include <string>
#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>

#include "include/solver.h"
#include "operator.h"
#include "var.h"

using Eigen::SparseMatrix;
using Eigen::VectorXd;

class Latent {
protected:
    unsigned n_obs;
    VectorXd Theta, Mu, W, Grad;

    SparseMatrix<double,0,int> A;
    Operator K;
    Var var;

public:
    Latent() {}
    // Latent(Rcpp::List){};
    ~Latent() {}

    void init_var(Rcpp::List var_in) {
        string type       = Rcpp::as<string>     (var_in["type"]);
        Rcpp::List v_init = Rcpp::as<Rcpp::List> (var_in["v_init"]);

        if (type == "ind_IG") {
            double a = Rcpp::as<double>  (v_init["a"]);
            double b = Rcpp::as<double>  (v_init["b"]);
            var = ind_IG(n_obs, a, b);
        }
    }

    virtual void sample_V() {};

    // sample V given w and Y
    virtual void sample_cond_V() {};

    // compute the grad. wrt parameter
    virtual void compute_grad() {}; 

    void setTheta(VectorXd& theta) {Theta = theta; } 
    VectorXd& getTheta() {return Theta; } 
    VectorXd& getGrad() {return Grad; } 

    VectorXd& getMu() { return Mu; }
    VectorXd& getW()  { return W; }
    VectorXd& getV()  { return var.getV(); }

    SparseMatrix<double, 0, int>& getA()     { return A; }
    SparseMatrix<double, 0, int>& getK()     { return K.getK(); }
    SparseMatrix<double, 0, int>& get_dK()   { return K.get_dK(); }
    SparseMatrix<double, 0, int>& get_d2K()  { return K.get_d2K(); }
};

#endif