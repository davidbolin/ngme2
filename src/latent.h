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
    unsigned n_reg, n_paras; //regressors, parameters
    
    // data
    double mu, sigma;
    VectorXd W, h;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;

public:
    Latent(Rcpp::List latent_in) 
    : n_paras ( Rcpp::as< unsigned > (latent_in["n_paras"]) ), 
      n_reg   ( Rcpp::as< unsigned > (latent_in["n_reg"]) ),
      A       ( Rcpp::as< SparseMatrix<double,0,int> > (latent_in["A"])),
      mu      (1),
      sigma   (1),
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

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_reg; } 
    unsigned getThetaSize() const             {return n_paras; } 
    SparseMatrix<double, 0, int>& getA()      {return A; }
    
    const VectorXd& getW()  const             {return W; }
    void            setW(const VectorXd& W)   {this->W = W; }
    
    // Related to optimizer
    const VectorXd getTheta() const;
    const VectorXd getGrad();
    void           setTheta(const VectorXd);


    /*  2 Variance component   */
    const VectorXd& getV() const { return var->getV(); }
    virtual void sample_cond_V()=0;

    /*  3 Operator component   */
    SparseMatrix<double, 0, int>& getK()    { return ope->getK(); }
    SparseMatrix<double, 0, int>& get_dK()  { return ope->get_dK(); }
    SparseMatrix<double, 0, int>& get_d2K() { return ope->get_d2K(); }

    const double getKappa() const             {return ope->getKappa(); } 
    void         setKappa(const double kappa) {ope->setKappa(kappa);} 
    virtual double _grad_kappa()=0;
};

// according to the config, decide what parameters to return
inline const VectorXd Latent::getTheta() const {
    // for now, just return kappa
    VectorXd theta (1);
    theta << getKappa();
    return theta;
}

inline const VectorXd Latent::getGrad() {
    // for now, just return kappa
    VectorXd grad (1);
    grad << _grad_kappa();
    return grad;
}

inline void Latent::setTheta(const VectorXd v) {
    // for now, just set kappa
    setKappa(v(0));
}

#endif