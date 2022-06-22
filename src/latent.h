/*
    latent class:
        dim. of paramter: n_ope_params + 3
    // goal: change mu, sigma -> multi-var case
*/

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
    bool debug;
    int n_mesh, n_params, n_ope, n_var {1}; // n_params=n_ope + n_mu + n_sigma + n_var

    // indicate optimize (K, mu, sigma, var)
    int opt_flag[4] {1, 1, 1, 1};
    
    bool use_precond {false}, numer_grad {false}, use_num_hess {true};

    // mu and sigma
    MatrixXd B_mu,  B_sigma;
    VectorXd theta_mu, theta_sigma;
    
    int n_mu, n_sigma;
    VectorXd mu, sigma;

    // eps for numerical gradient.
    double trace, trace_eps, eps;

    VectorXd W, prevW, h;
    SparseMatrix<double,0,int> A;
    
    Operator *ope;
    Var *var;

    // solver
    lu_sparse_solver solver_K;
    cholesky_solver  solver_Q; // Q = KT diag(1/SV) K

public:
    Latent(Rcpp::List latent_in) 
    : debug         (Rcpp::as< bool >        (latent_in["debug"])),
      n_mesh        (Rcpp::as< int >         (latent_in["n_mesh"])),
    
      B_mu          (Rcpp::as< MatrixXd >    (latent_in["B_mu"])),
      B_sigma       (Rcpp::as< MatrixXd >    (latent_in["B_sigma"])),
      
      theta_mu      (Rcpp::as< VectorXd >    (latent_in["theta_mu"])),
      theta_sigma   (Rcpp::as< VectorXd >    (latent_in["theta_sigma"])),
      
      n_mu          (theta_mu.size()),
      n_sigma       (theta_sigma.size()),
    //   n_mu          (Rcpp::as< double >    (latent_in["n_mu"])),
    //   n_sigma       (Rcpp::as< double >    (latent_in["n_sigma"])),
      
      trace         (0),
      trace_eps     (0),
      eps           (0.001), 
      
      W             (n_mesh),
      prevW         (n_mesh),
      h             (Rcpp::as< VectorXd >                     (latent_in["h"])),
      A             (Rcpp::as< SparseMatrix<double,0,int> >   (latent_in["A"]))
    {
if (debug) std::cout << "constructor of latent" << std::endl;
        // general input
            string var_type = Rcpp::as<string>     (latent_in["var_type"]);
        
        Rcpp::List control_f = Rcpp::as<Rcpp::List> (latent_in["control_f"]);
            opt_flag[0]   = Rcpp::as<int>        (control_f["opt_operator"]);
            opt_flag[1]   = Rcpp::as<int>        (control_f["opt_mu"]);
            opt_flag[2]   = Rcpp::as<int>        (control_f["opt_sigma"]);
            opt_flag[3]   = Rcpp::as<int>        (control_f["opt_var"]);

            use_precond  = Rcpp::as<bool>        (control_f["use_precond"] );
            use_num_hess = Rcpp::as<bool>        (control_f["use_num_hess"]);
            numer_grad   = Rcpp::as<bool>        (control_f["numer_grad"]) ;
            eps          = Rcpp::as<double>      (control_f["eps"]) ;
            
        // init values
            set_theta_mu(theta_mu);
            set_theta_sigma(theta_sigma);
            
            // sigma        = Rcpp::as<double>      (control_f["init_sigma"]);

        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (latent_in["operator_in"]);
            n_ope = ope_in["n_params"];
            n_var = 1;
            n_params = n_ope + n_mu + n_sigma + n_var;

        // construct var
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (latent_in["var_in"]);
        if (var_type == "nig") {
            var = new ind_IG(var_in, n_mesh, h);
        } else if (var_type == "normal") {
            var = new normal(var_in, n_mesh, h);
            // Not optimizing mu
            opt_flag[1] = 0;  
        }
if (debug) std::cout << "finish constructor of latent" << std::endl;
    }
    ~Latent() {}

    /*  1 Model itself   */
    unsigned getSize() const                  {return n_mesh; } 
    unsigned get_n_params() const             {return n_params; } 
    SparseMatrix<double, 0, int>& getA()      {return A; }
    
    const VectorXd& getW()  const             {return W; }
    const VectorXd& getPrevW()  const         {return prevW; }
    void            setW(const VectorXd& W)   { prevW = this->W; this->W = W; }

    VectorXd getMean() const { 
        VectorXd mean;
        if (n_mu==1) {
            mean = mu(0) * (getV()-h);
        }
        else {
            mean = mu.cwiseProduct(getV()-h);
        }
        return mean; 
    }

    /*  2 Variance component   */
    const VectorXd& getV()     const { return var->getV(); }
    const VectorXd& getPrevV() const { return var->getPrevV(); }
    VectorXd getSV() const { 
        VectorXd V=getV(); 
        // if (n_sigma==1)
        //     return (V * pow(sigma(0),2)); 
        // else 
        return (sigma.cwiseProduct(sigma).cwiseProduct(V));
    }
    VectorXd getPrevSV() const { 
        VectorXd V=getPrevV(); 
        // if (n_sigma==1)
        //     return (V * pow(sigma(0),2)); 
        // else 
            return (sigma.cwiseProduct(sigma).cwiseProduct(V));
    }

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

    // virtual double get_theta_sigma() const        { return log(sigma); }
    // virtual void   set_theta_sigma(double theta)  { this->sigma = exp(theta); }

    // Parameter: var
    virtual double get_theta_var() const   { return var->get_theta_var(); }
    virtual void   set_theta_var(double v) { var->set_theta_var(v); }
    virtual double grad_theta_var()        { return var->grad_theta_var(); }

    // Output
    virtual Rcpp::List get_estimates() const=0;
};


/*    Optimizer related    */
inline const VectorXd Latent::get_parameter() const {
if (debug) std::cout << "Start latent get parameter"<< std::endl;   
    int n_ope = ope->get_n_params();
    
    VectorXd parameter (n_params);
        parameter.segment(0, n_ope)             = ope->get_parameter();
        parameter.segment(n_ope, n_mu)          = get_theta_mu();
        parameter.segment(n_ope+n_mu, n_sigma)  = get_theta_sigma();
        parameter(n_ope+n_mu+n_sigma)           = get_theta_var();
    
if (debug) std::cout << "Finish latent get parameter"<< std::endl;   
    return parameter;
}

inline const VectorXd Latent::get_grad() {
if (debug) std::cout << "Start latent gradient"<< std::endl;   
    int n_ope = ope->get_n_params();
    VectorXd grad (n_params);
auto grad1 = std::chrono::steady_clock::now();

    if (opt_flag[0]) grad.segment(0, n_ope)             = grad_K_parameter();     else grad.segment(0, n_ope) = VectorXd::Constant(n_ope, 0);
    if (opt_flag[1]) grad.segment(n_ope, n_mu)          = grad_theta_mu();        else grad.segment(n_ope, n_mu) = VectorXd::Constant(n_mu, 0);
    if (opt_flag[2]) grad.segment(n_ope+n_mu, n_sigma)  = grad_theta_sigma();     else grad.segment(n_ope+n_mu, n_sigma) = VectorXd::Constant(n_sigma, 0);
    if (opt_flag[3]) grad(n_ope+n_mu+n_sigma)           = grad_theta_var();       else grad(n_ope+n_mu+n_sigma) = 0;

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

    // if (opt_flag[0])  set_theta_kappa(theta(0)); 
    if (opt_flag[0])  set_K_parameter   (theta.segment(0, n_ope));
    if (opt_flag[1])  set_theta_mu      (theta.segment(n_ope, n_mu)); 
    if (opt_flag[2])  set_theta_sigma   (theta.segment(n_ope+n_mu, n_sigma)); 
    if (opt_flag[3])  set_theta_var     (theta(n_ope+n_mu+n_sigma)); 
}

// sigma>0 -> theta=log(sigma)
// return the gradient wrt. theta, theta=log(sigma)
inline VectorXd Latent::grad_theta_sigma() {
    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd prevV = getPrevV();

    VectorXd result(n_sigma);

    if (n_sigma < 0) {
        // stationary case
        // double msq = (K*W - mu(0)*(V-h)).cwiseProduct(V.cwiseInverse()).dot(K*W - mu(0)*(V-h));
        double msq = (K*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().dot(V.cwiseInverse());
        double grad = - n_mesh / sigma(0) + pow(sigma(0), -3) * msq;

        // hessian using prevous V
        // double msq2 = (K*prevW - mu(0)*(prevV-h)).cwiseProduct(prevV.cwiseInverse()).dot(K*prevW - mu(0)*(prevV-h));
        double msq2 = (K*prevW - mu.cwiseProduct(prevV-h)).array().pow(2).matrix().dot(prevV.cwiseInverse());
        double hess = n_mesh / pow(sigma(0), 2) - 3 * pow(sigma(0), -4) * msq2;
        
        // grad. wrt theta
        result(0) =  grad / (hess * sigma(0) + grad);
        result(0) = -1.0 / n_mesh * grad * sigma(0);
    } else {

        // double msq = (K*W - mu.cwiseProduct(V-h)).cwiseProduct(V.cwiseInverse()).dot(K*W - mu(0)*(V-h));
        // VectorXd vsq = (K*W - mu.cwiseProduct(V-h)).array().pow(2);
        VectorXd vsq = (K*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
        VectorXd grad (n_sigma);

        for (int l=0; l < n_sigma; l++) {
            VectorXd tmp1 = vsq.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Constant(n_mesh, 1);
            VectorXd tmp2 = B_sigma.col(l).cwiseProduct(tmp1);
            grad(l) = tmp2.sum();
        }

        MatrixXd hess (n_sigma, n_sigma);
        VectorXd tmp3 = -2*vsq.cwiseProduct(sigma.array().pow(-2).matrix());
        // change
        hess = B_sigma.transpose() * tmp3.asDiagonal() * B_sigma;

        result = - 1.0 / n_mesh * grad;
        // result = hess.llt().solve(grad);
    }

    return result;
}

inline VectorXd Latent::grad_theta_mu() {
if (debug) std::cout << "Start mu gradient"<< std::endl;   
    VectorXd result(n_mu);

    SparseMatrix<double> K = getK();
    VectorXd V = getV();
    VectorXd inv_V = V.cwiseInverse();
    VectorXd prevV = getPrevV();
    VectorXd prev_inv_V = prevV.cwiseInverse();
    
    if (n_mu == 1 && n_sigma == 1) {
        // stationary case
        double grad = pow(sigma(0),-2) * (V-h).cwiseProduct(inv_V).dot(K*W - mu(0) * (V-h));
        double hess = -pow(sigma(0),-2) * (prevV-h).cwiseProduct(prev_inv_V).dot(prevV-h);

        result(0) = grad / hess;
    }
    else {
        VectorXd grad (n_mu);
        for (int l=0; l < n_mu; l++) {
            grad(l) = (V-h).cwiseProduct( B_mu.col(l).cwiseQuotient(getSV()) ).dot(K*W - mu.cwiseProduct(V-h));
        }

        result = - 1.0 / n_mesh * grad;
    }

if (debug) {
// std::cout << "grad of mu=" << grad <<std::endl;
// std::cout << "hess of mu=" << hess <<std::endl;
}
    return result;
}

// only for stationary case, delete later
// W|V ~ N(K^-1 mu(V-h), sigma^2 K-1 diag(V) K-T)
inline double Latent::function_kappa(double eps) {
    SparseMatrix<double> K = ope->getK(0, eps);

    VectorXd V = getV();
    VectorXd SV = getSV();

    SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
    
    solver_Q.compute(Q);
    
    VectorXd tmp = K * W - mu.cwiseProduct(V-h);

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
    
    VectorXd tmp = K * W - mu.cwiseProduct(V-h);

    double l = 0.5 * solver_Q.logdet() 
               - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);

    return l;
}

// numerical gradient for K parameters
inline VectorXd Latent::numerical_grad() {
std::cout << "start numerical gradient" <<std::endl;
    int n_ope = ope->get_n_params();
    VectorXd params = ope->get_parameter();
    double val = function_K(params);

    VectorXd grad (n_ope);
    for (int i=0; i < n_ope; i++) {
        VectorXd params_add_eps = params;
            params_add_eps(i) += eps;
        double val_add_eps = function_K(params_add_eps);
        double num_g = (val_add_eps - val) / eps;
        
        if (!use_num_hess) {
            grad(i) = - num_g / n_mesh;
        } else {
            VectorXd params_minus_eps = params;
                params_minus_eps(i) -= eps;
            double val_minus_eps = function_K(params_minus_eps);
            double num_hess = (val_minus_eps + val_add_eps - 2*val) / pow(eps, 2);
            grad(i) = num_g / num_hess;
        }
    } 
    return grad;
}


#endif