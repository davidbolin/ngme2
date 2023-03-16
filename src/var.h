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
    bool fix_V, fix_nu, init_V;
public:
    Var() {}
    Var(const Rcpp::List& noise_list, unsigned long seed) :
        var_rng       (seed),
        noise_type    (Rcpp::as<string>  (noise_list["noise_type"])),
        nu            (Rcpp::as<double>  (noise_list["nu"])),
        h             (Rcpp::as<VectorXd>  (noise_list["h"])),
        n             (Rcpp::as<int>     (noise_list["n_noise"])),
        V             (n),
        prevV         (n),
        fix_V         (Rcpp::as<bool>    (noise_list["fix_V"])),
        fix_nu        (Rcpp::as<bool>    (noise_list["fix_nu"])),
        init_V        (Rcpp::as<bool>    (noise_list["init_V"]))
    {
// std::cout << "n = " << n << std::endl;
        if (noise_type == "normal" || !init_V) {
            V = VectorXd::Ones(n);
            prevV = VectorXd::Ones(n);
            fix_nu = true;
            fix_V = true;
        } else { // nig or gal
            if (!Rf_isNull(noise_list["V"])) {
// std::cout << "init V with 1" << std::endl;
                V = Rcpp::as< VectorXd > (noise_list["V"]);
                prevV = V;
            } else {
// std::cout << "sample V " << std::endl;
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
        if (noise_type == "normal") return;

        if (noise_type == "gal") {
            prevV = V;
            VectorXd nu_vec = VectorXd::Constant(n, nu);
            VectorXd zero_vec = VectorXd::Constant(n, 1e-14);
            if (!fix_V) V = rGIG_cpp(h * nu, 2 * nu_vec, zero_vec, var_rng());
        } else {  // nig and normal+nig
            prevV = V;
            VectorXd nu_vec = VectorXd::Constant(n, nu);
            if (!fix_V) V = rGIG_cpp(VectorXd::Constant(n, -0.5), nu_vec, nu_vec.cwiseProduct(h.cwiseProduct(h)), var_rng());

        }
    }

    void sample_cond_V(const VectorXd& a_inc_vec, const VectorXd& b_inc_vec) {
        if (noise_type == "normal") return;

        if (noise_type == "gal") {
            prevV = V;
            VectorXd p_vec = (nu * h) - VectorXd::Constant(n, 0.5);
            VectorXd a_vec = VectorXd::Constant(n, 2 * nu) + a_inc_vec;
            VectorXd b_vec = b_inc_vec;
            if (!fix_V) V = rGIG_cpp(p_vec, a_vec, b_vec, var_rng());
            // here a_inc_vec is given as (mu/sigma)^2 (latent.h)
            // sampling from GIG(hv-0.5, 2nu + a_inc_vec, 0)

            // if (fix_V) return;
            // prevV = V;
            // for (int i=0; i < n; i++) {
            //     V(i) = R::rgamma(h(i)*nu - 0.5, 1 / (nu + a_inc_vec(i)/2));
            // }
        } else { // nig and normal+nig
            prevV = V;
            VectorXd p_vec = VectorXd::Constant(n, -1);
            VectorXd a_vec = VectorXd::Constant(n, nu) + a_inc_vec;
            VectorXd b_vec = VectorXd::Constant(n, nu).cwiseProduct(h).cwiseProduct(h) + b_inc_vec;
            if (!fix_V) V = rGIG_cpp(p_vec, a_vec, b_vec, var_rng());
        }
    }

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

    double grad_log_nu() const {
        if (fix_nu || noise_type == "normal") return 0;

        // grad of log nu
        double grad = 0;
        if (noise_type == "gal") {
            for(int i=0; i < n; i++) {
                double nu_hi = nu * h[i];
                //digamma(0.1) = digamma(1.1) - 1/0.1;
                if(nu_hi > 1){
                    grad -=  nu_hi * R::digamma(nu_hi);
                } else {
                    grad -=  nu_hi * R::digamma(nu_hi + 1) - 1.;
                }
            grad += nu_hi * (1 - log(1/nu) + log(V(i))) - nu * V(i);
// std::cout << "grad in var = " << grad << std::endl;
           }
           grad = - grad / n;
        } else { // type == nig or normal+nig
            // df/dnu = 0.5 (2h + 1/nu - h^2/V - V)
            // df/d(log nu) = df/dnu * nu
            VectorXd tmp = 0.5 * (2*h + VectorXd::Constant(n, 1/nu)
                - h.cwiseProduct(h).cwiseQuotient(V) - V);
            double grad_nu = tmp.mean() * nu;

            VectorXd tmp2 = 0.5 * (2*h + VectorXd::Constant(n, 1/nu)
                - h.cwiseProduct(h).cwiseQuotient(prevV) - prevV);
            double grad_nu2 = tmp2.mean();

            grad = grad_nu * nu;
            double hess_nu = -0.5 * pow(nu, -2);
            // hess of log nu
            double hess = nu * grad_nu2 + nu * nu * hess_nu;

            // grad = -grad;        // not use hessian
            grad = grad / hess;     // use hessian
        }

        // if (abs(grad) > 2) { // bad V
            // return 0;
        // } else
        return grad;
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
