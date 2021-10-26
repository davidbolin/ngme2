#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include<eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Latent
{
private:
    VectorXd theta_K, theta_m, theta_V;
    VectorXd m, V;
    
    MatrixXd K; // function?
public:
    Latent(VectorXd K, VectorXd V) : K(K), V(V) {}

    VectorXd mean(); // K, m
    MatrixXd cov(); // V, K
};


#endif