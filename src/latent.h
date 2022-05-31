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
    const int latent_para {4};

    bool debug;
    int n_reg; //regressors
    
    // indicate optimize (kappa, mu, sigma, var)
    int opt_flag[4] {1, 1, 1, 1};
    
    bool use_precond {false}, numer_grad {false};

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
    : debug   ( Rcpp::as< bool >        (latent_in["debug"])),
      n_reg   ( Rcpp::as< unsigned >    (latent_in["n_reg"]) ),
      
      mu        (0),
      sigma     (1),
      trace     (0),
      trace_eps (0),
      eps       (0.001), 
      
      W       (n_reg),
      prevW   (n_reg),
      h       (Rcpp::as< VectorXd > (latent_in["h"])),
      A       (Rcpp::as< SparseMatrix<double,0,int> > (latent_in["A"]))
    {
        // Init opt. flag
        opt_flag[0]   = Rcpp::as<int>        (latent_in["opt_kappa"]);
        opt_flag[1]   = Rcpp::as<int>        (latent_in["opt_mu"]);
        opt_flag[2]   = Rcpp::as<int>        (latent_in["opt_sigma"]);
        opt_flag[3]   = Rcpp::as<int>        (latent_in["opt_var"]);

        use_precond = Rcpp::as<bool>        (latent_in["use_precond"]);
        numer_grad  = Rcpp::as<bool>        (latent_in["numer_grad"]);

        // Init var
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (latent_in["var_in"]);
        string type       = Rcpp::as<string>     (var_in["type"]);

        // Set initial values
        Rcpp::List init_value = Rcpp::as<Rcpp::List> (latent_in["init_value"]);
        mu           = Rcpp::as<double>  (init_value["mu"]);
        sigma        = Rcpp::as<double>  (init_value["sigma"]);
        
        // construct var
        if (type == "ind_IG") {
            double nu    = Rcpp::as<double>  (init_value["nu"]);
            var = new ind_IG(n_reg, nu, h);
        } else if (type == "normal") {
            var = new normal(n_reg, h);
            // Not optimizing mu
            opt_flag[1] = 0;  
        }
    }
    ~Latent() {}

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_reg; } 
    unsigned getThetaSize() const             {return latent_para; } 
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
        var->sample_cond_V(getK(), W, mu, sigma);
    }

    /*  3 Operator component   */
    SparseMatrix<double, 0, int>& getK()    { return ope->getK(); }
    SparseMatrix<double, 0, int>& get_dK()  { return ope->get_dK(); }
    SparseMatrix<double, 0, int>& get_d2K() { return ope->get_d2K(); }

    // Paramter kappa
    void   setKappa(double kappa) {
        ope->setKappa(kappa);
        
        if (!numer_grad) compute_trace(); 
    } 

    /* 4 for optimizer */
    const VectorXd getTheta() const;
    const VectorXd getGrad();
    void           setTheta(const VectorXd&);
    void           finishOpt(int i) {opt_flag[i] = 0; }

    // Parameter: kappa
    virtual double get_theta_kappa() const=0;
    virtual void   set_theta_kappa(double v)=0;
    virtual double grad_theta_kappa()=0;

    virtual double function_kappa(double eps);    
    
    void compute_trace() {
        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
// compute trace
        solver_K.computeKTK(K);

// auto timer_trace = std::chrono::steady_clock::now();
        SparseMatrix<double> M = dK;
        trace = solver_K.trace(M);
// std::cout << "time for the trace (ms): " << since(timer_trace).count() << std::endl;   

        // update trace_eps if using hessian
        if ((!numer_grad) && (use_precond)) {
            SparseMatrix<double> K = ope->getK(eps);
            SparseMatrix<double> dK = ope->get_dK(eps);
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
inline const VectorXd Latent::getTheta() const {
    VectorXd theta (latent_para);

    theta(0) = get_theta_kappa();
    theta(1) = get_mu();         
    theta(2) = get_theta_sigma();
    theta(3) = get_theta_var();  
    
    return theta;
}

inline const VectorXd Latent::getGrad() {
    VectorXd grad (latent_para);
auto grad1 = std::chrono::steady_clock::now();
    if (opt_flag[0]) grad(0) = grad_theta_kappa();         else grad(0) = 0;
    if (opt_flag[1]) grad(1) = grad_mu();                  else grad(1) = 0;
    if (opt_flag[2]) grad(2) = grad_theta_sigma();         else grad(2) = 0;
    if (opt_flag[3]) grad(3) = grad_theta_var();           else grad(3) = 0;

// DEBUG: checking grads
if (debug) {
    std::cout << "grad_kappa (ms): " << since(grad1).count() << std::endl;   
    std::cout << "******* grad of kappa is: " << grad(0) << std::endl;   
    std::cout << "******* grad of mu is:    " << grad(1) << std::endl;   
    std::cout << "******* grad of sigma is: " << grad(2) << std::endl;   
    std::cout << "******* grad of var   is: " << grad(3) << std::endl;
}
    return grad;
}

inline void Latent::setTheta(const VectorXd& theta) {
    if (opt_flag[0])  set_theta_kappa(theta(0)); 
    if (opt_flag[1])  set_mu(theta(1)); 
    if (opt_flag[2])  set_theta_sigma(theta(2)); 
    if (opt_flag[3])  set_theta_var(theta(3)); 
}

// sigma>0 -> theta=log(sigma)
// return the gradient wrt. theta, theta=log(sigma)
inline double Latent::grad_theta_sigma() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd prevV = getPrevV();

    double msq = (K*W - mu*(V-h)).cwiseProduct(V.cwiseInverse()).dot(K*W - mu*(V-h));
    double msq2 = (K*prevW - mu*(prevV-h)).cwiseProduct(prevV.cwiseInverse()).dot(K*prevW - mu*(prevV-h));

    double grad = - n_reg / sigma + pow(sigma, -3) * msq;

    // hessian using prevous V
    double hess = n_reg / pow(sigma, 2) - 3 * pow(sigma, -4) * msq2;
    
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

// W|V ~ N(K^-1 mu(V-h), sigma^2 K-1 diag(V) K-T)
inline double Latent::function_kappa(double eps) {
    SparseMatrix<double> K = ope->getK(eps);

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


#endif