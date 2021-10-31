#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include <eigen/Dense>
#include "operator.h"
#include "var.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Latent
{
private:
    VectorXd theta_K, theta_m, theta_V;

    VectorXd m;
    MatrixXd A; // comes from data ()

// 2 important components - Operator and Variance component
    Operator K;
    Var V;
public:
    // initFromList(Rlist)
    Latent(MatrixXd A) {}
    
    // sample V given w and Y
    virtual void sample_V() {}

    // VectorXd mean(); // K, m
    // MatrixXd cov(); // V, K
};


#endif