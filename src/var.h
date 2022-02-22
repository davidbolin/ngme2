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
    VectorXd V;
public: 
    Var(){}
    Var(Rcpp::List);
    ~Var(){}

    const VectorXd& getV() const {return V;}
    
    virtual void sample_V()=0;
    virtual void sample_cond_V(SparseMatrix<double>& K,
                                VectorXd& W,
                                VectorXd& h, 
                                double mu, 
                                double sigma)=0;
};


class ind_IG : public Var{
private:
    double a, b;
public:
    ind_IG(){}
    ind_IG(unsigned n, double a, double b) : a(a), b(b) {
        this->n = n;
        V.resize(n);
        sample_V();
    }

    void sample_V() {
        gig sampler;
        for(int i = 0; i < n; i++)
            V[i] = sampler.sample(-0.5, a, b);
    };

    // sample V given W
    // V|W, Y ~ GIG(p-0.5, a+mu, b+(K W + h mu)^2)
    void sample_cond_V(SparseMatrix<double>& K,
                       VectorXd& W,
                       VectorXd& h, 
                       double mu, 
                       double sigma
                       ) {

        // VectorXd arg_2 = VectorXd::Constant(n, a) + Mu.cwiseProduct(Mu)/(sigma*sigma);   //+ eta * VectorXd::Ones(temporal.rows());
        VectorXd p_vec = VectorXd::Constant(n, -1);
        VectorXd a_vec = VectorXd::Constant(n, a+std::pow((mu/sigma), 2));   //+ eta * VectorXd::Ones(temporal.rows());
        VectorXd b_vec = VectorXd::Constant(n, b) + std::pow(sigma, -2) * (K*W + mu*h).cwiseProduct((K * W + mu*h));
        V = rGIG_cpp(p_vec, a_vec, b_vec);

// std::cout << "cond V|W=" << V << std::endl;
    };
};

#endif