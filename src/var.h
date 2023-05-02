// nu = nu

#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include "sample_rGIG.h"

using Eigen::VectorXd;
using Eigen::SparseMatrix;
using std::string;

class Var {
private:
    std::mt19937 var_rng;

    string noise_type; // normal or nig
    double nu;
    VectorXd h;

    unsigned n;
    VectorXd V, prevV;
    bool fix_V, fix_nu, init_V, single_V, hessian;
public:
    Var() {}
    Var(const Rcpp::List& noise_list, unsigned long seed) :
        var_rng       (seed),
        noise_type    (Rcpp::as<string>   (noise_list["noise_type"])),
        nu            (Rcpp::as<double>   (noise_list["nu"])),
        h             (Rcpp::as<VectorXd> (noise_list["h"])),
        n             (Rcpp::as<int>      (noise_list["n_noise"])),
        V             (h),
        prevV         (h),
        fix_V         (Rcpp::as<bool>     (noise_list["fix_V"])),
        fix_nu        (Rcpp::as<bool>     (noise_list["fix_nu"])),
        init_V        (Rcpp::as<bool>     (noise_list["init_V"])),
        single_V      (Rcpp::as<bool>     (noise_list["single_V"])),
        hessian       (Rcpp::as<bool>     (noise_list["hessian"]))
    {
// std::cout << "n = " << n << std::endl;
        if (noise_type == "normal" || !init_V) {
            fix_nu = true;
            fix_V = true;
        } else { // nig or gal
            if (!Rf_isNull(noise_list["V"])) {
                V = Rcpp::as< VectorXd > (noise_list["V"]);
                prevV = V;
            } else {
                sample_V();
                sample_V();
            }
        }
    }
    ~Var() {}

    string get_noise_type() const {return noise_type;}
    const VectorXd& getV()     const {return V;}
    const VectorXd& getPrevV() const {return prevV;}

    void setPrevV(const VectorXd& V) { if (!fix_V) prevV = V; }

    void sample_V() {
        if (fix_V) return;

        if (noise_type == "normal") return;

        if (noise_type == "gal") {
            prevV = V;
            VectorXd nu_vec = VectorXd::Constant(n, nu);
            VectorXd zero_vec = VectorXd::Constant(n, 1e-14);
            V = rGIG_cpp(h * nu, 2 * nu_vec, zero_vec, var_rng());
        } else {  // nig and normal+nig
            prevV = V;
            VectorXd nu_vec = VectorXd::Constant(n, nu);
            V = rGIG_cpp(VectorXd::Constant(n, -0.5), nu_vec, nu_vec.cwiseProduct(h.cwiseProduct(h)), var_rng());
        }
    }

    // called only by random effects, only one single V
//     void sample_cond_V(double a_inc, double b_inc, int dim = 1) {
//         if (fix_V) return;
//         double single_V;
//         if (noise_type == "normal") return;

//         if (noise_type == "gal") {
//             prevV = V;
//             double p = nu - 0.5 * dim;
//             double a = 2*nu + a_inc;
//             double b = b_inc;
//             single_V = rGIG_cpp(p, a, b, var_rng());
//         } else { // nig and normal+nig
//             prevV = V;
//             double p = -0.5 - 0.5*dim;
//             double a = nu + a_inc;
//             double b = nu + b_inc;
// if (a <= 0 || b <= 0) {
// std::cout << "p, a, b = " << p << ", " << a << ", " << b << std::endl;
// throw;
// }
//             single_V = rGIG_cpp(p, a, b, var_rng());
//         }

//         V = VectorXd::Constant(n, single_V);
//     }

    void sample_cond_V(const VectorXd& a_inc_vec, const VectorXd& b_inc_vec, int dim = 1, bool same = false);
    double grad_log_nu() const ;

    double get_nu() const {
        return nu;
    }

    double get_log_nu() const {
        if (noise_type == "normal")
            return 0;
        else
            return log(nu);
    }
    void set_log_nu(double theta) {
        if (noise_type != "normal") nu = exp(theta);
        // else doing nothing
    }
};
#endif

// a=b=nu
// class ind_IG : public Var{
// private:
//     double nu;
// public:
//     ind_IG() {}

//     ind_IG(double nu, unsigned n, unsigned long seed)
//     : nu(0)
//     {
//         this->seed = seed;
//         var_rng.seed(seed);
//         nu = nu;
//         this->n = n;

//         V.resize(n); prevV.resize(n);
//         sample_V(); sample_V(); // sample twice
//     }

//     double get_nu() const {return nu;}

//     // optimizer related
//     double get_log_nu() const   { return log(nu); }
//     void   set_log_nu(double theta) { nu = exp(theta); }
//     double grad_log_nu() const {
//         VectorXd tmp = VectorXd::Constant(n, 1+1/(2*nu)) - 0.5*V - VectorXd::Constant(n, 1).cwiseQuotient(2*V);
//         double grad = tmp.mean();
//         double hess = -0.5 * pow(nu, -2);

//         grad = grad / (hess * nu + grad);
//         return grad;
//     }

//     // sampling realted
//     void sample_V() {
//         prevV = V;

//         // gig sampler;
//         // for(int i = 0; i < n; i++)
//         //     V[i] = sampler.sample(-0.5, nu, nu);
//         VectorXd nu_vec = VectorXd::Constant(n, nu);
//         if (!fix_V) V = rGIG_cpp(VectorXd::Constant(n, -0.5), nu_vec, nu_vec, var_rng());
//     };

//     // sample cond. V
//     // V|W, Y ~ GIG(p-0.5, a + (mu/sigma)^2, b + ((KW + h mu) / sigma)^2) for process
//     // V|X, Y ~ GIG(p-0.5, a + (mu/sigma)^2, b + ((Y-Xb-AW) / sigma)^2) for block
//     // V|X, Y ~ GIG(p-0.5, a + a_inc_vec, b + b_inc_vec) for general use
//     void sample_cond_V(const VectorXd& a_inc_vec,
//                        const VectorXd& b_inc_vec
//     ) {
//         prevV = V;
//         // VectorXd arg_2 = VectorXd::Constant(n, a) + Mu.cwiseProduct(Mu)/(sigma*sigma);   //+ eta * VectorXd::Ones(temporal.rows());
//         // VectorXd a_vec = VectorXd::Constant(n, nu+std::pow((mu / sigma), 2));   //+ eta * VectorXd::Ones(temporal.rows());
//         // VectorXd b_vec = VectorXd::Constant(n, nu) + std::pow(sigma, -2) * (K*W + mu * h).cwiseProduct((K * W + mu * h));
//         // VectorXd b_vec = VectorXd::Constant(n, nu).array() + sigma.array().pow(-2).cwiseProduct( (K*W + mu.cwiseProduct(h)).array().pow(2) );

//         VectorXd p_vec = VectorXd::Constant(n, -1);
//         VectorXd a_vec = VectorXd::Constant(n, nu) + a_inc_vec;
//         VectorXd b_vec = VectorXd::Constant(n, nu) + b_inc_vec;
//         if (!fix_V) V = rGIG_cpp(p_vec, a_vec, b_vec, var_rng());
//     };
// };

// the case for normal noise
// can i get rid of h here??
// class normal : public Var {
// public:
//     normal() {}

//     normal(unsigned n) {
//         this->n = n;

//         V.resize(n); prevV.resize(n);
//         sample_V(); sample_V(); // sample twice
//     }

//     // normal(unsigned n, VectorXd h) {
//     //     this->n = n;
//     //     this->h = h;

//     //     V.resize(n); prevV.resize(n);
//     //     sample_V(); sample_V(); // sample twice
//     // }

//     double get_nu() const {return 0;}

//     // nothing to optimize
//     double get_log_nu() const   { return 1; }
//     void   set_log_nu(double theta) {}
//     double grad_log_nu() const {
//         return 0;
//     }

//     // return V=h
//     void sample_V() {
//         prevV = VectorXd::Ones(n);
//         V = VectorXd::Ones(n);
//     };

//     // return V=h
//     void sample_cond_V(const VectorXd& a_inc_vec,
//                        const VectorXd& b_inc_vec
//     ) {
//         prevV = V;
//         V = VectorXd::Ones(n);
//     };
// };
