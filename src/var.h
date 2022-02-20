#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "sample_rGIG.h"

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
    
    virtual void sample_V()=0;
    virtual void sample_cond_V(VectorXd&, 
                       VectorXd&, 
                       SparseMatrix<double>&)=0;
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
    void sample_cond_V(VectorXd& W, VectorXd& Mu, 
                       SparseMatrix<double>& K
                       ) {
        
        VectorXd h = VectorXd::Constant(n_obs, 1);
        
        double sigma = 1;
        VectorXd arg_1;
        arg_1.setOnes(n_obs);
        arg_1 = -arg_1;

        VectorXd arg_2 = VectorXd::Constant(n_obs, a) + Mu.cwiseProduct(Mu)/(sigma*sigma);   //+ eta * VectorXd::Ones(temporal.rows());
        VectorXd arg_3 = VectorXd::Constant(n_obs, b) + (K * W + h).cwiseProduct((K * W + h));
        V = rGIG_cpp(arg_1, arg_2, arg_3);
std::cout << "cond V|W=" << V << std::endl;

// std::cout << "arg_2" << arg_2 << std::endl;
// std::cout << "arg_3" << arg_3 << std::endl;
    };
};

#endif