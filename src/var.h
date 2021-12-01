#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <string>
#include "rgig.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::string;

class Var {
protected:
    unsigned n_obs;
    VectorXd V;
    VectorXd grad;
    
    MatrixXd hess;

public: 
    Var(){}
    Var(Rcpp::List);
    ~Var(){}

    VectorXd& const getV()  {return V;}
    
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

    // sample V given w and Y
    void sample_cond_V() {

    };
};

#endif
