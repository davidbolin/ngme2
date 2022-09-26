#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <random>
#include <string>
#include <cmath>
#include "sample_rGIG.h"

using Eigen::VectorXd;
using Eigen::SparseMatrix;
using std::string;

class Var {
protected:
    unsigned long seed;
    unsigned n;
    VectorXd V, prevV;
    std::mt19937 var_rng;
    bool fix_V {false};
public:
    Var() : V(n), prevV(n) {}

    const VectorXd& getV()     const {return V;}
    const VectorXd& getPrevV() const {return prevV;}

    void fixV() {
        fix_V = true;
    }

    void setV(const VectorXd& newV) {
        prevV = V;
        V = newV;
    }

    virtual void sample_V()=0;
    virtual void sample_cond_V(const VectorXd& a_inc_vec, const VectorXd& b_inc_vec)=0;

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
    ind_IG() {}

    ind_IG(double theta_V, unsigned n, unsigned long seed)
    : nu(0)
    {
        this->seed = seed;
        var_rng.seed(seed);
        nu = theta_V;
        this->n = n;

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

        // gig sampler;
        // for(int i = 0; i < n; i++)
        //     V[i] = sampler.sample(-0.5, nu, nu);
        VectorXd nu_vec = VectorXd::Constant(n, nu);
        if (!fix_V) V = rGIG_cpp(VectorXd::Constant(n, -0.5), nu_vec, nu_vec, var_rng());
    };

    // sample cond. V
    // V|W, Y ~ GIG(p-0.5, a + (mu/sigma)^2, b + ((KW + h mu) / sigma)^2) for process
    // V|X, Y ~ GIG(p-0.5, a + (mu/sigma)^2, b + ((Y-Xb-AW) / sigma)^2) for block
    // V|X, Y ~ GIG(p-0.5, a + a_inc_vec, b + b_inc_vec) for general use
    void sample_cond_V(const VectorXd& a_inc_vec,
                       const VectorXd& b_inc_vec
    ) {
        prevV = V;
        // VectorXd arg_2 = VectorXd::Constant(n, a) + Mu.cwiseProduct(Mu)/(sigma*sigma);   //+ eta * VectorXd::Ones(temporal.rows());
        // VectorXd a_vec = VectorXd::Constant(n, nu+std::pow((mu / sigma), 2));   //+ eta * VectorXd::Ones(temporal.rows());
        // VectorXd b_vec = VectorXd::Constant(n, nu) + std::pow(sigma, -2) * (K*W + mu * h).cwiseProduct((K * W + mu * h));
        // VectorXd b_vec = VectorXd::Constant(n, nu).array() + sigma.array().pow(-2).cwiseProduct( (K*W + mu.cwiseProduct(h)).array().pow(2) );

        VectorXd p_vec = VectorXd::Constant(n, -1);
        VectorXd a_vec = VectorXd::Constant(n, nu) + a_inc_vec;
        VectorXd b_vec = VectorXd::Constant(n, nu) + b_inc_vec;
        if (!fix_V) V = rGIG_cpp(p_vec, a_vec, b_vec, var_rng());
    };
};

// the case for normal noise
// can i get rid of h here??
class normal : public Var {
public:
    normal() {}

    normal(unsigned n) {
        this->n = n;

        V.resize(n); prevV.resize(n);
        sample_V(); sample_V(); // sample twice
    }

    // normal(unsigned n, VectorXd h) {
    //     this->n = n;
    //     this->h = h;

    //     V.resize(n); prevV.resize(n);
    //     sample_V(); sample_V(); // sample twice
    // }

    double get_var() const {return 0;}

    // nothing to optimize
    double get_theta_var() const   { return 1; }
    void   set_theta_var(double theta) {}
    double grad_theta_var() const {
        return 0;
    }

    // return V=h
    void sample_V() {
        prevV = VectorXd::Ones(n);
        V = VectorXd::Ones(n);
    };

    // return V=h
    void sample_cond_V(const VectorXd& a_inc_vec,
                       const VectorXd& b_inc_vec
    ) {
        prevV = V;
        V = VectorXd::Ones(n);
    };
};

#endif