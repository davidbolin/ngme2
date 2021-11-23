#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include <string>
#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include "operator.h"
#include "var.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Latent {
// friend class Var;
protected:
    unsigned n_obs;
    VectorXd theta_K, theta_m, theta_V;

    // VectorXd m = 0;
    MatrixXd A; // comes from data ()

    Operator K;
    Var var;
public:
    Latent() {}
    // Latent(Rcpp::List){};
    ~Latent() {}

    void init_var(Rcpp::List var_in) {
        string type = var_in["type"];
        Rcpp::List v_init = var_in["v_init"];

        if (type == "ind_IG") {
            double a = v_init["a"];
            double b = v_init["b"];
            var = ind_IG(n_obs, a, b);
        }
    }

    virtual void sample_V() {};

    // sample V given w and Y
    virtual void sample_cond_V() {};


    MatrixXd getA() {
        return A;
    }

    VectorXd getV() {
        return var.getV();
    }

    MatrixXd getK()  {
        return K.getK();
    }
};

class AR : public Latent {
// class AR  {
// friend class Var;

public:
    AR(){}
    AR(Rcpp::List ar1_in) {
        A = ar1_in["A"];
        n_obs = ar1_in["n"];

        Rcpp::List ope_in = ar1_in["operator_in"]; // containing C and G
        K = GC(ope_in);

        Rcpp::List var_in = ar1_in["var_in"];
        init_var(var_in);
    }

    void sample_V() {
    }

    void sample_cond_V() {};
};

#endif