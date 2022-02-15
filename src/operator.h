#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

using Eigen::SparseMatrix;
using Eigen::VectorXd;

class Operator {
protected:
    VectorXd theta_K;
    SparseMatrix<double, 0, int> K, dK, d2K;
public:
    Operator() {};
    Operator(VectorXd theta_K) {};

    virtual void update(VectorXd& theta_K)=0;

    // getter for K, dK, d2K
    SparseMatrix<double, 0, int>& getK()    {return K;}
    SparseMatrix<double, 0, int>& get_dK()  {return dK;}
    SparseMatrix<double, 0, int>& get_d2K() {return d2K;}
};

// fit for AR and Matern 

class GC : public Operator {
private:
    double theta {1};
    SparseMatrix<double, 0, int> G, C;

public:
    GC(Rcpp::List ope_in) {
        theta = ope_in["a_init"];
        G = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]);
        C = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]);
        
        K = theta * G + C;
        dK = C;
        d2K = 0 * C;
    }

    void update(VectorXd& theta_K) {
        assert (theta_K.size() == 1);
        K = theta_K(0) * G + C;
    }

};

#endif