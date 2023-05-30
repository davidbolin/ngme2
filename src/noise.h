// nu = nu

#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include "sample_rGIG.h"

using Eigen::VectorXd;
using std::string;

class NoiseUtil {
public:
    // update gig representation (p, a, b)
    static void update_gig(
        const string& noise_type,
        double nu,
        Eigen::Ref<Eigen::VectorXd> p,
        Eigen::Ref<Eigen::VectorXd> a,
        Eigen::Ref<Eigen::VectorXd> b,
        const VectorXd& h
    );

    static double grad_theta_nu(
        const string& noise_type,
        double nu,
        const VectorXd& V,
        const VectorXd& prevV,
        const VectorXd& h
    );

    // wrapper of sample GIG with extra argument
    static void sample_V(
        Eigen::Ref<Eigen::VectorXd> V,
        const string& noise_type,
        const VectorXd& p,
        const VectorXd& a,
        const VectorXd& b,
        std::mt19937& rng,
        bool single_V = false
    ) {
        if (single_V) {
            V = VectorXd::Constant(1, rGIG_cpp(p[0], a[0], b[0], rng()));
        } else {
            V = rGIG_cpp(p, a, b, rng());
        }
    }

    // ---- h = 1 ----
    static void update_gig(
        const string& noise_type,
        double nu,
        Eigen::Ref<Eigen::VectorXd> p,
        Eigen::Ref<Eigen::VectorXd> a,
        Eigen::Ref<Eigen::VectorXd> b
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

    static VectorXd rnorm_vec(int n, double mu, double sigma, unsigned long seed=0) {
        std::mt19937 norm_rng(seed);
        std::normal_distribution<double> rnorm {0,1};
        VectorXd out(n);
        for (int i = 0; i < n; i++)
        {
            // out[i] = R::rnorm(mu, sigma);
            out[i] = rnorm(norm_rng) * sigma + mu;
        }
        return (out);
    }
};

#endif
