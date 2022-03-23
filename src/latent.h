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
      A       ( Rcpp::as< SparseMatrix<double,0,int> > (latent_in["A"])),
      mu      (0),
      sigma   (1),
      W       (n_reg),
      prevW   (n_reg),
      h       (VectorXd::Constant(n_reg, 1))
    {
        // Init Data
        
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
    virtual double grad_theta_var()        { return var->grad_theta_var();}

    // sigma
    virtual double get_theta_sigma() const        { return log(sigma); }
    virtual void   set_theta_sigma(double theta)  { this->sigma = exp(theta); }
    virtual double grad_theta_sigma();

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
    if (opt_kappa) grad(0) = grad_theta_kappa(); else grad(0) = 0;
std::cout << "grad_kappa (ms): " << since(grad1).count() << std::endl;   

    if (opt_mu)    grad(1) = grad_mu();          else grad(1) = 0;

    if (opt_sigma) grad(2) = grad_theta_sigma(); else grad(2) = 0;

    if (opt_var)   grad(3) = grad_theta_var();   else grad(3) = 0;

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

// compute trace in 1 way
auto timer_init = std::chrono::steady_clock::now();
    solver_K.computeKKT(K);
std::cout << "time for compute K (ms): " << since(timer_init).count() << std::endl;

auto timer_trace = std::chrono::steady_clock::now();
    SparseMatrix<double> M = dK.transpose() * K;
    double trace = solver_K.trace(M);
std::cout << "time for the trace (ms): " << since(timer_trace).count() << std::endl;   
std::cout << "trace1=: " << trace << std::endl;   

// 2. ordinary way
// auto timer_init2 = std::chrono::steady_clock::now();
//     SparseLU<SparseMatrix<double>> eigen_solver_K;
//     eigen_solver_K.compute(K);
// std::cout << "time for compute K2 (ms): " << since(timer_init2).count() << std::endl;

// ///// takes 1s to compute the trace
// auto timer_trace2 = std::chrono::steady_clock::now();
//     SparseMatrix<double> trace_mat = eigen_solver_K.solve(dK);
    
//     double trace2 = trace_mat.diagonal().sum();
// std::cout << "time for the trace2 (ms): " << since(timer_trace2).count() << std::endl;   
// std::cout << "trace2=: " << trace2 << std::endl;   
    
    // tr(dK * K^-1)
    // double lhs = solver_K.trace(dK); 
    // double rhs = W.transpose() * dK.transpose() * 
    //             (VectorXd::Constant(n_reg, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * mu);
    double rhs = (dK*W).cwiseProduct(SV.cwiseInverse()).cwiseProduct(K * W + (h - V) * mu).sum();
    return (rhs - trace) / n_reg;
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

inline double Latent::grad_mu() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd inv_V = V.cwiseInverse();
    VectorXd Vmh = V-h;
    
    VectorXd prevV = getPrevV();
    VectorXd prev_inv_V = prevV.cwiseInverse();
    VectorXd prevVmh = prevV-h;

    // double hess_mu = -(Vmh).transpose() * inv_SV.asDiagonal() * Vmh;  // get previous V
    // double g = (Vmh).transpose() * inv_SV.asDiagonal() * (K*W - mu*Vmh);
    double hess_mu = -pow(sigma,-2) * (prevVmh).cwiseProduct(prev_inv_V).cwiseProduct(prevVmh).sum();
    double g = pow(sigma,-2) * (Vmh).cwiseProduct(inv_V).cwiseProduct(K*W-mu*Vmh).sum();

    return g / hess_mu;
}

// inline double Latent::grad_mu() {
//     SparseMatrix<double> K = getK();
//     solver_K.compute(K);
//     VectorXd mean = getMean();
//     VectorXd m = solver_K.solve(mean);

//     VectorXd V = getV();
//         VectorXd inv_V = VectorXd::Constant(V.size(), 1).cwiseQuotient(V);
//     VectorXd SV = getSV();
//         VectorXd inv_SV = VectorXd::Constant(SV.size(), 1).cwiseQuotient(SV);

//     SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
//     SparseMatrix<double> QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;
//         if (!Q_analyzed) { solver_Q.analyze(QQ); Q_analyzed=true; }
//     solver_Q.compute(QQ);


//     VectorXd tmp = K.transpose() * inv_SV.asDiagonal() * getMean() + 
//                    pow(sigma_eps, -2) * A.transpose() * Y;
    
//     VectorXd mm = solver_Q.solve(tmp);
//         tmp = V - h;
//         VectorXd tmp2 = A * solver_K.solve(tmp);
//     double g = tmp2.transpose() * inv_V.asDiagonal() * (Y - mu*tmp2 - A * (mm-m));

//     return -g/n_reg;
// }

#endif