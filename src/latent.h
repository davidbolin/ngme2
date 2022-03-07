#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include <string>
#include <iostream>
#include <cmath>
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
    int n_reg, n_paras {4}; //regressors, parameters
    
    // indicate which parameter to optimize
    bool opt_mu {false}, opt_sigma {false}, opt_kappa {false}, opt_var {false};

    // data
    VectorXd Y;
    double sigma_eps;

    double mu, sigma;
    VectorXd W, h;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;
    lu_sparse_solver solver_K;

public:
    Latent(Rcpp::List latent_in) 
    : n_reg   ( Rcpp::as< unsigned > (latent_in["n_reg"]) ),
      A       ( Rcpp::as< SparseMatrix<double,0,int> > (latent_in["A"])),
      mu      (1),
      sigma   (1),
      W       (n_reg),
      h       (VectorXd::Constant(n_reg, 1)),
      sigma_eps(1)
    {
        // Init opt. flag
        opt_mu = Rcpp::as<bool>         (latent_in["opt_mu"]);
        opt_sigma = Rcpp::as<bool>      (latent_in["opt_sigma"]);
        opt_kappa = Rcpp::as<bool>      (latent_in["opt_kappa"]);
        opt_var = Rcpp::as<bool>        (latent_in["opt_var"]);

        // Init var
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (latent_in["var_in"]);
        string type       = Rcpp::as<string>     (var_in["type"]);
        Rcpp::List v_init = Rcpp::as<Rcpp::List> (var_in["v_init"]);

        if (type == "ind_IG") {
            double nu = Rcpp::as<double>  (v_init["nu"]);
            var = new ind_IG(n_reg, nu);
        }
    }
    ~Latent() {}

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_reg; } 
    unsigned getThetaSize() const             {return n_paras; } 
    SparseMatrix<double, 0, int>& getA()      {return A; }
    
    const VectorXd& getW()  const             {return W; }
    void            setW(const VectorXd& W)   {this->W = W; }

    const VectorXd getMean() const { return mu * (getV() - h); }

    // Parameter mu
    double getMu() const     {return mu;} 
    void   setMu (double mu) {this->mu = mu;} 

    // virtual double get_theta_kappa() const=0;
    // virtual void   set_theta_kappa(const VectorXd& v)=0;
    // virtual double _grad_theta_kappa()=0;

    /*  2 Variance component   */
    VectorXd getSV() const { VectorXd V=getV(); return (V*sigma); } // local object
    const VectorXd& getV() const { return var->getV(); }
    virtual void sample_cond_V()=0;

    /*  3 Operator component   */
    SparseMatrix<double, 0, int>& getK()    { return ope->getK(); }
    SparseMatrix<double, 0, int>& get_dK()  { return ope->get_dK(); }
    SparseMatrix<double, 0, int>& get_d2K() { return ope->get_d2K(); }

    // Paramter kappa
    double getKappa() const       {return ope->getKappa(); } 
    void   setKappa(double kappa) {ope->setKappa(kappa);} 
    
    /* 4 for optimizer */

    const VectorXd getTheta() const;
    const VectorXd getGrad();
    void           setTheta(const VectorXd&);

    // kappa
    virtual double get_theta_kappa() const=0;
    virtual void   set_theta_kappa(double v)=0;
    virtual double grad_theta_kappa()=0;
    double _grad_kappa();
    
    // nu
    virtual double get_theta_var() const   { return var->get_theta_var(); }
    virtual void   set_theta_var(double v) { var->set_theta_var(v); }
    virtual double grad_theta_var()        { return var->grad_theta_var();}

    // sigma
    virtual double get_theta_sigma() const        { return log(sigma); }
    virtual void   set_theta_sigma(double theta)  { this->sigma = exp(theta); }
    virtual double grad_theta_sigma();

};

/*    Optimizer related    */
//  bool opt_mu {false}, opt_sigma {false}, opt_kappa {false}, opt_var {false};
inline const VectorXd Latent::getTheta() const {
    VectorXd theta (n_paras);
    if (opt_kappa) theta(0) = get_theta_kappa(); else theta(0) = 0;
    if (opt_mu)    theta(1) = getMu();           else theta(1) = 0;
    if (opt_sigma) theta(2) = get_theta_sigma(); else theta(2) = 0;
    if (opt_var)   theta(3) = get_theta_var();   else theta(3) = 0;
    
    return theta;
}

inline const VectorXd Latent::getGrad() {
    VectorXd grad (n_paras);
    if (opt_kappa) grad(0) = grad_theta_kappa(); else grad(0) = 0;
    if (opt_mu)    grad(1) = 0;                  else grad(1) = 0;
    if (opt_sigma) grad(2) = grad_theta_sigma(); else grad(2) = 0;
    if (opt_var)   grad(3) = grad_theta_var();   else grad(3) = 0;

    return grad;
}

inline void Latent::setTheta(const VectorXd& theta) {
    if (opt_kappa) set_theta_kappa(theta(0)); 
    // if (opt_sigma) set_theta_sigma(theta(2)); 
    if (opt_sigma) set_theta_sigma(theta(2)); 
    if (opt_var)   set_theta_sigma(theta(3)); 
    
}

inline double Latent::_grad_kappa() {
    SparseMatrix<double> K = getK();
    SparseMatrix<double> dK = get_dK();
    VectorXd V = getV();
    solver_K.compute(K);

    double lhs = solver_K.trace(dK); // tr(dK * K^-1)
    // 2. Compute the rest
    double rhs = W.transpose() * dK.transpose() * 
                (VectorXd::Constant(n_reg, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * mu);
    return (rhs - lhs) / n_reg;
}

// sigma>0 -> theta=log(sigma)
// return the gradient wrt. theta, theta=log(sigma)
inline double Latent::grad_theta_sigma() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
        VectorXd inv_V = VectorXd::Constant(V.size(), 1).cwiseQuotient(V);

    VectorXd mm = (K*W + (h-V)*mu);

    double l = log(sigma) * (-n_reg);
    double r = (mm.transpose() * inv_V.asDiagonal() * mm); r = r * pow(sigma, -3);

    // grad. wrt theta
    return -(l+r)/n_reg * sigma;
}

#endif