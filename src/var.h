#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <cmath>
#include "sample_rGIG.h"

using Eigen::VectorXd;
using Eigen::SparseMatrix;
using std::string;

class Var {
protected:
    unsigned n;
    VectorXd V, prevV, h;
public: 
    Var() : V(n), prevV(n) {}
    ~Var(){}

    const VectorXd& getV()     const {return V;}
    const VectorXd& getPrevV() const {return prevV;}
    
    virtual void sample_V()=0;
    virtual void sample_cond_V(SparseMatrix<double>& K,
                                VectorXd& W,
                                VectorXd& mu, 
                                double sigma)=0;

    virtual double get_var()               const=0;
    virtual double get_theta_var()         const=0;
    virtual double grad_theta_var()        const=0;
    virtual void   set_theta_var(double)=0;
};

// a=b=nu
class ind_IG : public Var{
private:
    double nu;
public:
    ind_IG(){}
    
    ind_IG(Rcpp::List var_in, unsigned n, VectorXd h) 
    : nu(0)
    {   
        nu = (Rcpp::as<double>) (var_in["nu"]);
        this->n = n; 
        this->h = h;

        V.resize(n); prevV.resize(n);
        sample_V(); sample_V(); // sample twice
    }

    double get_var() const {return nu;} 
    
    // optimizer related
    double get_theta_var() const   { return log(nu); }
    void   set_theta_var(double theta) { nu = exp(theta); }
    
    double grad_theta_var() const {
        VectorXd tmp = VectorXd::Constant(n, 1+1/(2*nu)) - 0.5*V - VectorXd::Constant(n, 1).cwiseQuotient(2*V);
        double grad = tmp.mean();
        double hess = -0.5 * pow(nu, -2);
        
        grad = grad / (hess * nu + grad);
        return grad;
    }

    // sampling realted
    void sample_V() {
        prevV = V;

        gig sampler;
        for(int i = 0; i < n; i++)
            V[i] = sampler.sample(-0.5, nu, nu);
    };

    // sample V given W
    // V|W, Y ~ GIG(p-0.5, a+mu, b+(K W + h mu)^2)
    void sample_cond_V(SparseMatrix<double>& K,
                       VectorXd& W,
                       VectorXd& mu, 
                       double sigma
                       ) {
        prevV = V;
        if (mu.size() == 1) {
        // VectorXd arg_2 = VectorXd::Constant(n, a) + Mu.cwiseProduct(Mu)/(sigma*sigma);   //+ eta * VectorXd::Ones(temporal.rows());
            VectorXd p_vec = VectorXd::Constant(n, -1);
            VectorXd a_vec = VectorXd::Constant(n, nu+std::pow((mu(0) / sigma), 2));   //+ eta * VectorXd::Ones(temporal.rows());
            VectorXd b_vec = VectorXd::Constant(n, nu) + std::pow(sigma, -2) * (K*W + mu(0) * h).cwiseProduct((K * W + mu(0) * h));
            V = rGIG_cpp(p_vec, a_vec, b_vec);
        } else {
            throw ("sample conditional V not implemented for non-stationary case");
        }
    };
};

// adding the case for normal
class normal : public Var {
public:
    normal() {}
    normal(Rcpp::List var_in, unsigned n, VectorXd h) {
        this->n = n;
        this->h = h;
        
        V.resize(n); prevV.resize(n);
        sample_V(); sample_V(); // sample twice
    }

    double get_var() const {return 0;} 

    // nothing to optimize
    double get_theta_var() const   { return 0; }
    void   set_theta_var(double theta) {}    
    double grad_theta_var() const {
        return 0;
    }

    // return V=h
    void sample_V() {
        prevV = h;
        V = h;
    };

    // return V=h
    void sample_cond_V(SparseMatrix<double>& K,
                       VectorXd& W,
                       VectorXd& mu, 
                       double sigma
                       ) {
        prevV = h;
        V = h;
    };
};

#endif