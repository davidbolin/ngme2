#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include <string>
#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>

#include "include/timer.h"
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

    double mu, sigma;
    VectorXd W, prevW, h;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;

    // solver
    lu_sparse_solver solver_K;
    cholesky_solver  solver_Q;
    bool Q_analyzed {false};

public:
    Latent(Rcpp::List latent_in) 
    : n_reg   ( Rcpp::as< unsigned > (latent_in["n_reg"]) ),
      
      mu      (0),
      sigma   (1),
      
      W       (n_reg),
      prevW   (n_reg),
      h       (VectorXd::Constant(n_reg, 1)),
      A       (Rcpp::as< SparseMatrix<double,0,int> > (latent_in["A"]))
    {
        // Init opt. flag
        opt_mu = Rcpp::as<bool>         (latent_in["opt_mu"]);
        opt_sigma = Rcpp::as<bool>      (latent_in["opt_sigma"]);
        opt_kappa = Rcpp::as<bool>      (latent_in["opt_kappa"]);
        opt_var = Rcpp::as<bool>        (latent_in["opt_var"]);

        // Init var
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (latent_in["var_in"]);
        string type       = Rcpp::as<string>     (var_in["type"]);
        // Set initial values
        Rcpp::List init_value = Rcpp::as<Rcpp::List> (latent_in["init_value"]);
        mu           = Rcpp::as<double>  (init_value["mu"]);
        sigma        = Rcpp::as<double>  (init_value["sigma"]);
        double nu    = Rcpp::as<double>  (init_value["nu"]);
        if (type == "ind_IG") {
            var = new ind_IG(n_reg, nu);
        }

    }
    ~Latent() {}

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_reg; } 
    unsigned getThetaSize() const             {return n_paras; } 
    SparseMatrix<double, 0, int>& getA()      {return A; }
    
    const VectorXd& getW()  const             {return W; }
    void            setW(const VectorXd& W)   { prevW = this->W; this->W = W; }

    VectorXd getMean() const { return mu * (getV() - h); }

    /*  2 Variance component   */
    VectorXd getSV() const { VectorXd V=getV(); return (V*pow(sigma,2)); }
    const VectorXd& getV()     const { return var->getV(); }
    const VectorXd& getPrevV() const { return var->getPrevV(); }
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
    virtual double grad_theta_var()        { 
        
        // sigma -> sigma.tilde
        // double nu = var->get_var();
        // double dst_dnu = -0.5 * (sigma * mu * mu) / pow(1+nu*mu*mu, 1.5);
        double dst_dnu = 1;
        
        return var->grad_theta_var() * dst_dnu;
    }

    /* Old parameterization */
    virtual double get_theta_sigma() const        { return log(sigma); }
    virtual void   set_theta_sigma(double theta)  { this->sigma = exp(theta); }
    virtual double grad_theta_sigma();

    /* Use new sigma tilde parametrization */
    // sigma.tilde = sigma / sqrt(1 + eta * nu^2)
    // double get_sigma_tilde() const {return sigma / sqrt(1 + mu*mu * var->get_var()); }
    // virtual double get_theta_sigma() const        { return log(get_sigma_tilde()); }
    // virtual void   set_theta_sigma(double theta)  { this->sigma = exp(theta) * sqrt(1 + mu*mu * var->get_var()); }
    // virtual double grad_theta_sigma();

    // mu
    double get_mu() const     {return mu;} 
    void   set_mu (double mu) {this->mu = mu;} 
    virtual double grad_mu();
};

/*    Optimizer related    */
inline const VectorXd Latent::getTheta() const {
    VectorXd theta (n_paras);

    theta(0) = get_theta_kappa();
    theta(1) = get_mu();         
    theta(2) = get_theta_sigma();
    theta(3) = get_theta_var();  
    
    return theta;
}

inline const VectorXd Latent::getGrad() {
    VectorXd grad (n_paras);
auto grad1 = std::chrono::steady_clock::now();
    if (opt_kappa) grad(0) = grad_theta_kappa() / n_reg; else grad(0) = 0;

    if (opt_mu)    grad(1) = grad_mu();                  else grad(1) = 0;

    if (opt_sigma) grad(2) = grad_theta_sigma() / n_reg; else grad(2) = 0;

    if (opt_var)   grad(3) = grad_theta_var();           else grad(3) = 0;

std::cout << "grad_kappa (ms): " << since(grad1).count() << std::endl;   
    return grad;
}

inline void Latent::setTheta(const VectorXd& theta) {
    if (opt_kappa) set_theta_kappa(theta(0)); 
    if (opt_mu)    set_mu(theta(1)); 
    if (opt_sigma) set_theta_sigma(theta(2)); 
    if (opt_var)   set_theta_var(theta(3)); 
}

inline double Latent::_grad_kappa() {
    SparseMatrix<double> K = getK();
    SparseMatrix<double> dK = get_dK();
    VectorXd V = getV();
    VectorXd SV = getSV();

// compute trace
auto timer_init = std::chrono::steady_clock::now();
    solver_K.computeKKT(K);
std::cout << "time for compute K (ms): " << since(timer_init).count() << std::endl;

auto timer_trace = std::chrono::steady_clock::now();
    SparseMatrix<double> M = dK.transpose() * K;
    double trace = solver_K.trace(M);
std::cout << "time for the trace (ms): " << since(timer_trace).count() << std::endl;   
// std::cout << "trace1=: " << trace << std::endl;   

double rhs = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V) * mu);
    return (rhs - trace);
}

// sigma>0 -> theta=log(sigma)
// return the gradient wrt. theta, theta=log(sigma)
inline double Latent::grad_theta_sigma() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd prevV = getPrevV();

    double msq = (K*W - mu*(V-h)).cwiseProduct(V.cwiseInverse()).dot(K*W - mu*(V-h));
    double msq2 = (K*W - mu*(prevV-h)).cwiseProduct(prevV.cwiseInverse()).dot(K*W - mu*(prevV-h));

    double grad = - n_reg / sigma + pow(sigma, -3) * msq;

    // hessian using prevous V
    double hess = n_reg / pow(sigma, 2) - 3 * pow(sigma, -4) * msq2;
    
    // grad. wrt theta
std::cout << "******* grad of sigma is: " << grad / (hess * sigma + grad) * n_reg << std::endl;   

    return grad / (hess * sigma + grad) * n_reg;
}

inline double Latent::grad_mu() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd inv_V = V.cwiseInverse();
    
    VectorXd prevV = getPrevV();
    VectorXd prev_inv_V = prevV.cwiseInverse();

    // double hess_mu = -(Vmh).transpose() * inv_SV.asDiagonal() * Vmh;  // get previous V
    // double g = (Vmh).transpose() * inv_SV.asDiagonal() * (K*W - mu*Vmh);
    double hess = -pow(sigma,-2) * (prevV-h).cwiseProduct(prev_inv_V).dot(prevV-h);
    double grad = pow(sigma,-2) * (V-h).cwiseProduct(inv_V).dot(K*W - mu*(V-h));

    return grad / hess;
}

#endif