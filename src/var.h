#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "rgig.h"

using Eigen::VectorXd;
using Eigen::SparseMatrix;
using std::string;

class Var {
protected:
    unsigned n_obs;
    VectorXd V;
    VectorXd grad;
    
    SparseMatrix<double> hess;
public: 
    Var(){}
    Var(Rcpp::List);
    ~Var(){}

    VectorXd&  getV() {return V;}
    
    void sample_V() {};
    // sample V given w and Y
    void sample_cond_V() {};
};


class ind_IG : public Var{
private:
    double a, b;
public:
    ind_IG(){}
    ind_IG(unsigned n, double a, double b) : a(a), b(b) {
        n_obs = n;
        V.resize(n_obs);
        sample_V();
    }

    void sample_V() {
        gig sampler;
        for(int i = 0; i < n_obs; i++)
            V[i] = sampler.sample(-0.5, a, b);
    };

    // sample V given W
    // V|W ~ GIG(p-0.5, a+mu, b+(K W + h mu)^2)
    void sample_cond_V(VectorXd W) {
        VectorXd Mu;

        VectorXd arg_1;
        arg_1.setOnes(n_obs);
        arg_1 = -arg_1;

        VectorXd arg_2 = Mu.cwiseProduct(Mu)/(sigma*sigma) + eta * VectorXd::Ones(temporal.rows());
        VectorXd arg_3 = b + (K * W + h * Mu).cwiseProduct((K * W + h * Mu));

        // V = rGIG_cpp(arg_1, arg_2, arg_3);
    };
};

#endif