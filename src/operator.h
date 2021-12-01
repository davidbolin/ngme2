#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <cassert>

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class Operator {
protected:
    VectorXd theta_K;
    SparseMatrix<double, 0, int> K, dK, d2K;
public:
    Operator() {};
    Operator(VectorXd theta_K) {};

    virtual void update(VectorXd theta_K) {}; 

    // getter for K, dK, d2K
    SparseMatrix<double, 0, int>& const getK() {return K;}
    SparseMatrix<double, 0, int>& const get_dK() {return dK;}
    SparseMatrix<double, 0, int>& const get_d2K() {return d2K;}
};

// fit for AR and Matern 

class GC : public Operator {
private:
    double theta {1};
    MatrixXd G, C;

public:
    GC(Rcpp::List ope_in) {
        theta = ope_in["a_init"];
        G = ope_in["G"];
        C = ope_in["C"];
        
        K = G + theta * C;
        dK = C;
        d2K = 0 * C;
    }

    // GC(double theta) {
    //     K = G + theta * C;
    //     dK = C;
    //     d2K = 0 * C;
    // }

    void update(VectorXd theta_K) {
        assert (theta_K.size() == 1);
        
        K = G + theta_K(1) * C;
    }

};

#endif