// nu = nu

#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>

using Eigen::VectorXd;
using std::string;

class V_related {
public:
    // update gig representation (p, a, b)
    static void update_gig(
        const string& noise_type,
        double nu,
        VectorXd& p,
        VectorXd& a,
        VectorXd& b,
        const VectorXd& h
    );

    static double grad_theta_nu(
        const string& noise_type,
        double nu,
        const VectorXd& V,
        const VectorXd& prevV,
        const VectorXd& h
    );

    // h=1 version
    static void update_gig(
        const string& noise_type,
        double nu,
        VectorXd& p,
        VectorXd& a,
        VectorXd& b
    ) {
        VectorXd h = VectorXd::Ones(p.size());
        update_gig(noise_type, nu, p, a, b, h);
    }

    static double grad_theta_nu(
        const string& noise_type,
        double nu,
        const VectorXd& V,
        const VectorXd& prevV
    ) {
        VectorXd h = VectorXd::Ones(V.size());
        return grad_theta_nu(noise_type, nu, V, prevV, h);
    }
};

#endif
