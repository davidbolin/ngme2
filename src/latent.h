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
    int n_params;

    bool debug;
    int n_mesh; 

    // indicate optimize (kappa, mu, sigma, var)
    int opt_flag[4] {1, 1, 1, 1};
    
    bool use_precond {false}, numer_grad {false};

    // eps for numerical gradient.
    double mu, sigma, trace, trace_eps, eps;

    VectorXd W, prevW, h;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;

    // solver
    lu_sparse_solver solver_K;
    cholesky_solver  solver_Q; // Q = KT diag(1/SV) K

public:
    Latent(Rcpp::List latent_in) 
    : debug    ( Rcpp::as< bool >        (latent_in["debug"])),
      
      mu        (0),
      sigma     (1),
      trace     (0),
      trace_eps (0),
      eps       (0.001), 
      
      W       (n_mesh),
      prevW   (n_mesh),
      h       ( Rcpp::as< VectorXd >                     (latent_in["h"])),
      A       ( Rcpp::as< SparseMatrix<double,0,int> >   (latent_in["A"]))
    {
        // general input
            n_mesh   = Rcpp::as< int >                          (latent_in["n_mesh"]);
            n_params = Rcpp::as< int >                          (latent_in["n_la_params"]);

            string var_type = Rcpp::as<string>     (latent_in["var_type"]);
        
        Rcpp::List control_f = Rcpp::as<Rcpp::List> (latent_in["control_f"]);
            opt_flag[0]   = Rcpp::as<int>        (control_f["opt_operator"]);
            opt_flag[1]   = Rcpp::as<int>        (control_f["opt_mu"]);
            opt_flag[2]   = Rcpp::as<int>        (control_f["opt_sigma"]);
            opt_flag[3]   = Rcpp::as<int>        (control_f["opt_var"]);

            use_precond = Rcpp::as<bool>        (control_f["use_precond"]);
            numer_grad  = Rcpp::as<bool>        (control_f["numer_grad"]);
            eps         = Rcpp::as<double>      (control_f["eps"]);
            
            // init values
            mu           = Rcpp::as<double>     (control_f["init_mu"]);
            sigma        = Rcpp::as<double>     (control_f["init_sigma"]);
        
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (latent_in["var_in"]);
        
        // construct var
        if (var_type == "nig") {
            var = new ind_IG(var_in, n_mesh, h);
        } else if (var_type == "normal") {
            var = new normal(var_in, n_mesh, h);
            // Not optimizing mu
            opt_flag[1] = 0;  
        }
    }
    ~Latent() {}

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_mesh; } 
    unsigned get_n_params() const             {return n_params; } 
    SparseMatrix<double, 0, int>& getA()      {return A; }
    
    const VectorXd& getW()  const             {return W; }
    const VectorXd& getPrevW()  const         {return prevW; }
    void            setW(const VectorXd& W)   { prevW = this->W; this->W = W; }

    VectorXd getMean() const { return mu * (getV() - h); }

    /*  2 Variance component   */
    VectorXd getSV() const { VectorXd V=getV(); return (V*pow(sigma,2)); }
    const VectorXd& getV()     const { return var->getV(); }
    const VectorXd& getPrevV() const { return var->getPrevV(); }
    
    void sample_cond_V() {
        var->sample_cond_V(ope->getK(), W, mu, sigma);
    }

    /*  3 Operator component   */
    // for block use
    SparseMatrix<double, 0, int>& getK()    { return ope->getK(); }
    SparseMatrix<double, 0, int>& get_dK()  { return ope->get_dK(); }
    SparseMatrix<double, 0, int>& get_d2K() { return ope->get_d2K(); }

    /* 4 for optimizer */
    const VectorXd get_parameter() const;
    const VectorXd get_grad();
    void           set_parameter(const VectorXd&);
    void           finishOpt(int i) {opt_flag[i] = 0; }

    // Parameter: Operator
    virtual VectorXd    get_K_parameter() const {return ope->get_parameter(); } // no change of variable
    virtual VectorXd    grad_K_parameter() { return numerical_grad(); }
    virtual void        set_K_parameter(VectorXd params) {ope->set_parameter(params); }


    // used for ar
    virtual double function_kappa(double eps);    
    
    // used for general case
    virtual double function_K(VectorXd parameter);    
    virtual VectorXd numerical_grad(); // given eps

    // only for stationary case
    void compute_trace() {
        SparseMatrix<double> K = ope->getK();
        SparseMatrix<double> dK = ope->get_dK();
// compute trace
        solver_K.computeKTK(K);

// auto timer_trace = std::chrono::steady_clock::now();
        SparseMatrix<double> M = dK;
        trace = solver_K.trace(M);
// std::cout << "time for the trace (ms): " << since(timer_trace).count() << std::endl;   

        // update trace_eps if using hessian
        if ((!numer_grad) && (use_precond)) {
            SparseMatrix<double> K = ope->getK(0, eps);
            SparseMatrix<double> dK = ope->get_dK(0, eps);
            SparseMatrix<double> M = dK;

            solver_K.computeKTK(K);
            trace_eps = solver_K.trace(M);
        }
    };
    
    // Parameter: nu
    virtual double get_theta_var() const   { return var->get_theta_var(); }
    virtual void   set_theta_var(double v) { var->set_theta_var(v); }
    virtual double grad_theta_var()        { 
        return var->grad_theta_var();
    }

    // Parameter: sigma
    virtual double get_theta_sigma() const        { return log(sigma); }
    virtual void   set_theta_sigma(double theta)  { this->sigma = exp(theta); }
    virtual double grad_theta_sigma();

    // Parameter: mu
    double get_mu() const     {return mu;} 
    void   set_mu(double mu) {this->mu = mu;} 
    virtual double grad_mu();

    // Output
    virtual Rcpp::List get_estimates() const=0;
};


/*    Optimizer related    */
inline const VectorXd Latent::get_parameter() const {
    int n_ope = ope->get_n_params();
    
    VectorXd parameter (n_params);
        parameter.segment(0, n_ope) = ope->get_parameter();
        parameter(n_ope)            = get_mu();         
        parameter(n_ope+1)          = get_theta_sigma();
        parameter(n_ope+2)          = get_theta_var();  
    
    return parameter;
}

inline const VectorXd Latent::get_grad() {
    int n_ope = ope->get_n_params();
    VectorXd grad (n_params);
auto grad1 = std::chrono::steady_clock::now();

    if (opt_flag[0]) grad.segment(0, n_ope) = grad_K_parameter();     else grad.segment(0, n_ope) = VectorXd::Constant(n_ope, 0);
    if (opt_flag[1]) grad(n_ope)            = grad_mu();              else grad(n_ope) = 0;
    if (opt_flag[2]) grad(n_ope+1)          = grad_theta_sigma();     else grad(n_ope+1) = 0;
    if (opt_flag[3]) grad(n_ope+2)          = grad_theta_var();       else grad(n_ope+2) = 0;

// DEBUG: checking grads
if (debug) {
    // std::cout << "grad_kappa (ms): " << since(grad1).count() << std::endl;   
    // std::cout << "******* grad of kappa is: " << grad(0) << std::endl;   
    // std::cout << "******* grad of mu is:    " << grad(1) << std::endl;   
    // std::cout << "******* grad of sigma is: " << grad(2) << std::endl;   
    // std::cout << "******* grad of var   is: " << grad(3) << std::endl;
}
    return grad;
}

inline void Latent::set_parameter(const VectorXd& theta) {
    int n_ope = ope->get_n_params();

    // if (opt_flag[0])  set_theta_kappa(theta(0)); 
    if (opt_flag[0])  set_K_parameter(theta.segment(0, n_ope));
    if (opt_flag[1])  set_mu(theta(n_ope)); 
    if (opt_flag[2])  set_theta_sigma(theta(n_ope+1)); 
    if (opt_flag[3])  set_theta_var(theta(n_ope+2)); 
}

// sigma>0 -> theta=log(sigma)
// return the gradient wrt. theta, theta=log(sigma)
inline double Latent::grad_theta_sigma() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd prevV = getPrevV();

    double msq = (K*W - mu*(V-h)).cwiseProduct(V.cwiseInverse()).dot(K*W - mu*(V-h));
    double msq2 = (K*prevW - mu*(prevV-h)).cwiseProduct(prevV.cwiseInverse()).dot(K*prevW - mu*(prevV-h));

    double grad = - n_mesh / sigma + pow(sigma, -3) * msq;

    // hessian using prevous V
    double hess = n_mesh / pow(sigma, 2) - 3 * pow(sigma, -4) * msq2;
    
    // grad. wrt theta

    return grad / (hess * sigma + grad);
}


inline double Latent::grad_mu() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd inv_V = V.cwiseInverse();
    
    VectorXd prevV = getPrevV();
    VectorXd prev_inv_V = prevV.cwiseInverse();

    double grad = pow(sigma,-2) * (V-h).cwiseProduct(inv_V).dot(K*W - mu*(V-h));
    double hess = -pow(sigma,-2) * (prevV-h).cwiseProduct(prev_inv_V).dot(prevV-h);

    return grad / hess;
}

// only for stationary case, delete later
// W|V ~ N(K^-1 mu(V-h), sigma^2 K-1 diag(V) K-T)
inline double Latent::function_kappa(double eps) {
    SparseMatrix<double> K = ope->getK(0, eps);

    VectorXd V = getV();
    VectorXd SV = getSV();

    SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
    
    solver_Q.compute(Q);
    
    VectorXd tmp = K * W - mu*(V-h);

    double l = 0.5 * solver_Q.logdet() 
               - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
                // - 0.5 * (prevW-mean).transpose() * Q * (prevW-mean);

    return l;
}

// function_K(params += ( 0,0,eps,0,0) )
inline double Latent::function_K(VectorXd parameter) {
    assert(parameter.size()==ope->get_n_params());
    SparseMatrix<double> K = ope->getK(parameter);
    
    VectorXd V = getV();
    VectorXd SV = getSV();

    SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
    
    solver_Q.compute(Q);
    
    VectorXd tmp = K * W - mu*(V-h);

    double l = 0.5 * solver_Q.logdet() 
               - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);

    return l;
}

inline VectorXd Latent::numerical_grad() {
    int n_ope = ope->get_n_params();
    VectorXd params = ope->get_parameter();
    double val = function_K(params);

    VectorXd grad (n_ope);
    for (int i=0; i < n_ope; i++) {
        VectorXd params_add_eps = params;
            params_add_eps(i) += eps;
        double val_add_eps = function_K(params_add_eps);
        
        double num_g = (val_add_eps - val) / eps;
        
        if (!use_precond) {
            grad(i) = - num_g / n_mesh;
        } else {
            VectorXd params_minus_eps = params;
                params_minus_eps(i) -= eps;
            double val_minus_eps = function_K(params_minus_eps);
            double num_hess = (val_minus_eps + val_add_eps - 2*val) / pow(eps, 2);
            grad(i) = num_g / num_hess;
        }
    } 
}


#endif