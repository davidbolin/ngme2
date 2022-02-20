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
    unsigned n_reg, n_paras;
    double mu;

    VectorXd W, h, Theta;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;

public:
    // Latent() {}
    Latent(Rcpp::List latent_in) 
    : n_paras ( Rcpp::as< unsigned > (latent_in["n_paras"]) ), 
      n_reg   ( Rcpp::as< unsigned > (latent_in["n_reg"]) ),
      A       ( Rcpp::as< SparseMatrix<double,0,int> > (latent_in["A"])),
      mu      (1),
      Theta   (n_paras),
      W       (n_reg),
      h       (VectorXd::Constant(n_reg, 1))
    {
        // Init var
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (latent_in["var_in"]);
        string type       = Rcpp::as<string>     (var_in["type"]);
        Rcpp::List v_init = Rcpp::as<Rcpp::List> (var_in["v_init"]);

        if (type == "ind_IG") {
            double a = Rcpp::as<double>  (v_init["a"]);
            double b = Rcpp::as<double>  (v_init["b"]);
            var = new ind_IG(n_reg, a, b);
        }
    }
    ~Latent() {}

    virtual void sample_V()=0;
    virtual void sample_cond_V()=0;
    virtual VectorXd& getTheta() {return Theta; } 

    unsigned getThetaSize() {return n_paras; } 

    VectorXd& getW()  { return W; }
    VectorXd& getV()  { return var->getV(); }
    
    virtual void setTheta(VectorXd& theta) {Theta = theta; } 
    virtual void setW(VectorXd& new_W)  { W = new_W; }
    
    SparseMatrix<double, 0, int>& getA()     { return A; }
    SparseMatrix<double, 0, int>& getK()     { return ope->getK(); }
    SparseMatrix<double, 0, int>& get_dK()   { return ope->get_dK(); }
    SparseMatrix<double, 0, int>& get_d2K()  { return ope->get_d2K(); }
};

#endif