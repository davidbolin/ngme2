/*
    latent class:
    basically, providing get, grad, set (theta_K, ...)
        1. theta_K
        2. theta_mu
        3. theta_sigma
        4. theta_V
*/

#ifndef NGME_LATANT_H
#define NGME_LATANT_H

// #include<Eigen/IterativeLinearSolvers>
#include <string>
#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <random>

#include "include/timer.h"
#include "include/solver.h"
#include "operator.h"
#include "var.h"

using Eigen::SparseMatrix;
using Eigen::VectorXd;

enum Latent_fix_flag {
    latent_fix_ope, latent_fix_mu, latent_fix_sigma, 
    latent_fix_var, latent_fix_V, latent_fix_W
};
const int LATENT_FIX_FLAG_SIZE = 6;

class Latent {
protected:
    unsigned long seed;
    string model_type, noise_type;
    bool debug;
    int n_mesh, n_params, n_ope, n_var {1}; // n_params=n_ope + n_theta_mu + n_theta_sigma + n_var

    // indicate fixing (K, mu, sigma, var)  (1 means fix)
    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0};
    
    bool use_precond {false}, numer_grad {false};
    bool symmetricK {false};

    // mu and sigma
    MatrixXd B_mu,  B_sigma;
    VectorXd theta_mu, theta_sigma;
    
    int n_theta_mu, n_theta_sigma;
    VectorXd mu, sigma;

    // eps for numerical gradient.
    double trace, trace_eps, eps;

    VectorXd W, prevW, h;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;

    // solver
    cholesky_solver chol_solver_K;
    lu_sparse_solver lu_solver_K;

    bool use_iter_solver {false};
    iterative_solver CG_solver_K;

    cholesky_solver solver_Q; // Q = KT diag(1/SV) K
    
    std::mt19937 latent_rng;
public:
    Latent(Rcpp::List, unsigned long seed);
    ~Latent() {}

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_mesh; } 
    unsigned get_n_params() const             {return n_params; } 
    SparseMatrix<double, 0, int>& getA()      {return A; }
    
    const VectorXd& getW()  const             {return W; }
    const VectorXd& getPrevW()  const         {return prevW; }
    void            setW(const VectorXd& W)   { 
        if (!fix_flag[latent_fix_W]) {
            prevW = this->W;
            this->W = W; 
        }
    }

    VectorXd getMean() const { return mu.cwiseProduct(getV()-h); }
    VectorXd getMeanKX() { 
        SparseMatrix<double> K = getK();
        VectorXd res (n_mesh);
        if (!symmetricK) {
            lu_solver_K.compute(K);
            res = lu_solver_K.solve(getMean());
        } else {
            chol_solver_K.compute(K);
            res = chol_solver_K.solve(getMean());
        }
        return A * res;
    }

    /*  2 Variance component   */
    const VectorXd& getV()     const { return var->getV(); }
    const VectorXd& getPrevV() const { return var->getPrevV(); }
    VectorXd getSV() const {
        VectorXd V=getV(); 
        return (sigma.array().pow(2).matrix().cwiseProduct(V));
    }
    VectorXd getPrevSV() const { 
        VectorXd prevV=getPrevV(); 
        return (sigma.array().pow(2).matrix().cwiseProduct(prevV));
    }

    void sample_cond_V() {
        VectorXd tmp = (ope->getK() * W + mu.cwiseProduct(h));
        VectorXd a_inc_vec = mu.cwiseQuotient(sigma).array().pow(2);
        VectorXd b_inc_vec = tmp.cwiseQuotient(sigma).array().pow(2);
        var->sample_cond_V(a_inc_vec, b_inc_vec);
    }

    /*  3 Operator component   */
    // for block use
    SparseMatrix<double, 0, int>& getK()    { return ope->getK(); }
    // SparseMatrix<double, 0, int>& get_dK()  { return ope->get_dK(); }
    // SparseMatrix<double, 0, int>& get_d2K() { return ope->get_d2K(); }

    /* 4 for optimizer */
    const VectorXd get_parameter() const;
    const VectorXd get_grad();
    void           set_parameter(const VectorXd&);
    void           finishOpt(int i) {fix_flag[i] = 0; }

    // Parameter: Operator (override this if do change of variable)
    virtual VectorXd    get_theta_K() const {
        return ope->get_parameter(); 
    } 
    virtual VectorXd    grad_theta_K() { return numerical_grad(); }
    virtual void        set_theta_K(VectorXd params) {ope->set_parameter(params); }

    // deprecated
    // virtual double function_kappa(double eps);    
    virtual double function_K(VectorXd parameter);   
    
    // used for general case
    virtual double function_K(SparseMatrix<double> K);
    virtual VectorXd numerical_grad(); // given eps

    // update the trace value
    void compute_trace() {
        SparseMatrix<double> K = ope->getK();
        SparseMatrix<double> dK = ope->get_dK(0);
// compute trace
// auto timer_trace = std::chrono::steady_clock::now();

        SparseMatrix<double> M = dK;
        if (!use_iter_solver) {
            if (!symmetricK) {
                lu_solver_K.computeKTK(K);
                trace = lu_solver_K.trace(M);
            } else {
                chol_solver_K.compute(K);
                trace = chol_solver_K.trace(M);
            }
        } else {
            if (!symmetricK) {
                // BiCG solver
                throw("Not implemented yet");
            } else {
                SparseMatrix<double, RowMajor> KK = K;
                CG_solver_K.compute(K);
                trace = CG_solver_K.trace(M);
            }
        }

// std::cout << "trace ====== " << trace << std::endl;
// std::cout << "time for the trace (ms): " << since(timer_trace).count() << std::endl;   

        // update trace_eps if using hessian
        if ((!numer_grad) && (use_precond)) {
            SparseMatrix<double> K = ope->getK(0, eps);
            SparseMatrix<double> dK = ope->get_dK(0, 0, eps);
            SparseMatrix<double> M = dK;

        if (!use_iter_solver) {
            if (!symmetricK) {
                lu_solver_K.computeKTK(K);
                trace_eps = lu_solver_K.trace(M);
            } else {
                chol_solver_K.compute(K);
                trace_eps = chol_solver_K.trace(M);
            }
        } else {
            if (!symmetricK) {
                // BiCG solver
                throw("Not implemented yet");
            } else {
                CG_solver_K.compute(K);
                trace_eps = CG_solver_K.trace(M);
            }
        }
        }
    };

    // Parameter: mu
    VectorXd get_theta_mu() const {return theta_mu;} 
    void   set_theta_mu(VectorXd theta_mu)  {
        this->theta_mu = theta_mu;
        mu = (B_mu * theta_mu);
    } 
    virtual VectorXd grad_theta_mu();

    // Parameter: sigma
    virtual VectorXd get_theta_sigma() const { return theta_sigma; }
    virtual void set_theta_sigma(VectorXd theta_sigma) { 
        this->theta_sigma = theta_sigma;
        sigma = (B_sigma * theta_sigma).array().exp();
    }
    virtual VectorXd grad_theta_sigma();

    // Parameter: var
    virtual double get_theta_var() const   { return var->get_theta_var(); }
    virtual void   set_theta_var(double v) { var->set_theta_var(v); }
    virtual double grad_theta_var()        { return var->grad_theta_var(); }

    // Output
    virtual Rcpp::List get_estimates() const=0;
    
    Rcpp::List output() const {
        return Rcpp::List::create(
            Rcpp::Named("model_type")   = model_type,
            Rcpp::Named("noise_type")   = noise_type,
            Rcpp::Named("theta_mu")     = theta_mu,
            Rcpp::Named("theta_sigma")  = theta_sigma,
            Rcpp::Named("theta_V")      = var->get_theta_var(),
            Rcpp::Named("V")            = getV(),
            Rcpp::Named("W")            = W
        );
    }
};

/*    Optimizer related    */
inline const VectorXd Latent::get_parameter() const {
if (debug) std::cout << "Start latent get parameter"<< std::endl;   
    int n_ope = ope->get_n_params();
    
    VectorXd parameter (n_params);
        parameter.segment(0, n_ope)                         = get_theta_K();
        parameter.segment(n_ope, n_theta_mu)                = get_theta_mu();
        parameter.segment(n_ope+n_theta_mu, n_theta_sigma)  = get_theta_sigma();
        parameter(n_ope+n_theta_mu+n_theta_sigma)           = get_theta_var();
    
if (debug) std::cout << "End latent get parameter"<< std::endl;   
    return parameter;
}

inline const VectorXd Latent::get_grad() {
if (debug) std::cout << "Start latent gradient"<< std::endl;   
    int n_ope = ope->get_n_params();
    VectorXd grad (n_params);
    
auto grad1 = std::chrono::steady_clock::now();
    if (!fix_flag[latent_fix_ope])     grad.segment(0, n_ope)                        = grad_theta_K();         else grad.segment(0, n_ope) = VectorXd::Constant(n_ope, 0);
    if (!fix_flag[latent_fix_mu])      grad.segment(n_ope, n_theta_mu)               = grad_theta_mu();        else grad.segment(n_ope, n_theta_mu) = VectorXd::Constant(n_theta_mu, 0);
    if (!fix_flag[latent_fix_sigma])   grad.segment(n_ope+n_theta_mu, n_theta_sigma) = grad_theta_sigma();     else grad.segment(n_ope+n_theta_mu, n_theta_sigma) = VectorXd::Constant(n_theta_sigma, 0);
    if (!fix_flag[latent_fix_var])     grad(n_ope+n_theta_mu+n_theta_sigma)          = grad_theta_var();       else grad(n_ope+n_theta_mu+n_theta_sigma) = 0;

// DEBUG: checking grads
if (debug) {
    std::cout << "gradient time " << since(grad1).count() << std::endl;   
    std::cout << "gradient= " << grad << std::endl;   
}
    return grad;
}

inline void Latent::set_parameter(const VectorXd& theta) {
if (debug) std::cout << "Start latent set parameter"<< std::endl;   
    int n_ope = ope->get_n_params();

    set_theta_K       (theta.segment(0, n_ope));
    set_theta_mu      (theta.segment(n_ope, n_theta_mu)); 
    set_theta_sigma   (theta.segment(n_ope+n_theta_mu, n_theta_sigma)); 
    set_theta_var     (theta(n_ope+n_theta_mu+n_theta_sigma)); 
}

#endif