// nu = nu
#ifndef NGME_VAR_H
#define NGME_VAR_H

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include "sample_rGIG.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::string;

class NoiseUtil {
public:
    static constexpr double Pi = 3.14159265358979323846;
    // update gig representation (p, a, b)
    static void update_gig(
        const string& noise_type,
        const VectorXd& nu,
        Eigen::Ref<Eigen::VectorXd> p,
        Eigen::Ref<Eigen::VectorXd> a,
        Eigen::Ref<Eigen::VectorXd> b,
        const VectorXd& h,
        bool single_V = false
    );

    static VectorXd grad_theta_nu(
        const string& noise_type,
        const MatrixXd& B_nu,
        const VectorXd& nu,
        const VectorXd& V,
        const VectorXd& prevV,
        const VectorXd& h,
        double nu_lower_bound = 0.0,
        bool single_V = false
    );

    // Analytic Hessian wrt theta_nu for NIG (and normal_nig) case:
    // H = B_nu^T diag(nu ⊙ c) B_nu, where
    //   nu = nu_lower_bound + exp(B_nu theta_nu) (caller supplies effective 'nu'),
    //   c_i = h_i - h_i^2/(2 V_i) - V_i/2.
    // Returns negative Hessian (matching gradient sign convention).
    static MatrixXd hess_theta_nu(
        const string& noise_type,
        const MatrixXd& B_nu,
        const VectorXd& nu,
        const VectorXd& V,
        const VectorXd& h,
        double nu_lower_bound = 0.0,
        bool single_V = false
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

    // compute pi(V|theta_nu, h)
    static double log_density(
        const string& noise_type,
        const VectorXd& V,
        const VectorXd& h,
        const MatrixXd& B_nu,
        const VectorXd& theta_nu,
        bool single_V = false
    );

    static MatrixXd precond(
        double eps,
        const string& noise_type,
        const VectorXd& V,
        const VectorXd& h,
        const MatrixXd& B_nu,
        const VectorXd& theta_nu,
        bool single_V = false
    );

    // ---- h = 1 ----
    static void update_gig(
        const string& noise_type,
        const VectorXd& nu,
        Eigen::Ref<Eigen::VectorXd> p,
        Eigen::Ref<Eigen::VectorXd> a,
        Eigen::Ref<Eigen::VectorXd> b,
        bool single_V = false
    ) {
        VectorXd h = VectorXd::Ones(p.size());
        update_gig(noise_type, nu, p, a, b, h);
    }

    static VectorXd grad_theta_nu(
        const string& noise_type,
        const MatrixXd& B_nu,
        const VectorXd& nu,
        const VectorXd& V,
        const VectorXd& prevV,
        double nu_lower_bound = 0.0,
        bool single_V = false
    ) {
        VectorXd h = VectorXd::Ones(V.size());
        return grad_theta_nu(noise_type, B_nu, nu, V, prevV, h, nu_lower_bound, single_V);
    }

    // Convenience overload: Hessian with h = 1 vector
    static MatrixXd hess_theta_nu(
        const string& noise_type,
        const MatrixXd& B_nu,
        const VectorXd& nu,
        const VectorXd& V,
        double nu_lower_bound = 0.0,
        bool single_V = false
    ) {
        VectorXd h = VectorXd::Ones(V.size());
        return hess_theta_nu(noise_type, B_nu, nu, V, h, nu_lower_bound, single_V);
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
