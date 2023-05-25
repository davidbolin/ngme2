/*
    latent class:
    basically, providing get, grad, set (theta_K, ...)
        1. theta_K is unbounded
        2. theta_mu
        3. theta_sigma
        4. nu
*/

#ifndef NGME_LATANT_H
#define NGME_LATANT_H

// #include<Eigen/IterativeLinearSolvers>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <memory>

#include "include/timer.h"
#include "include/solver.h"
#include "operator.h"
#include "noise.h"

using std::exp;
using std::log;
using std::pow;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

enum Latent_fix_flag {
    latent_fix_theta_K,
    latent_fix_W,
    latent_fix_V,
    latent_fix_theta_mu,
    latent_fix_nu,
    latent_fix_theta_sigma
};
const int LATENT_FIX_FLAG_SIZE = 5;

class Latent {
protected:
    std::mt19937 latent_rng;
    string model_type;
    vector<string> noise_type;
    bool debug;
    int W_size, V_size, n_params, n_var {1}; // n_params=n_theta_K + n_theta_mu + n_theta_sigma + n_var

    // operator
    std::unique_ptr<Operator> ope, ope_add_eps;
    VectorXd h, theta_K;
    int n_theta_K;
    bool symmetricK, zero_trace, use_num_dK {false};

    vector<double> trace;
    double eps;

    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0};
    bool use_precond {false}, numer_grad {false};

    // mu and sigma, and sigma_normal (special case when using nig_normal case)
    MatrixXd B_mu, B_sigma, B_sigma_normal;
    VectorXd theta_mu, theta_sigma, theta_sigma_normal;

    // mu = Bmu * theta_mu
    // sigma = exp(Bsigma * theta_sigma)
    VectorXd mu, sigma, sigma_normal;
    int n_theta_mu, n_theta_sigma, n_nu, n_theta_sigma_normal;

    // for numerical gradient.
    VectorXd W, prevW, V, prevV;
    SparseMatrix<double,0,int> A;

    VectorXd p_vec, a_vec, b_vec, nu;
    bool single_V {false};

    // solver
    cholesky_solver chol_solver_K;
    lu_sparse_solver lu_solver_K;
    cholesky_solver solver_Q; // Q = KT diag(1/SV) K
public:
    Latent(const Rcpp::List&, unsigned long seed);
    ~Latent() {}

    void update_each_iter(bool init=false);

    /*  1 Model itself   */
    int get_W_size() const             {return W_size; }
    int get_V_size() const             {return V_size; }
    int get_n_params() const           {return n_params; }
    int get_n_theta_K() const          {return n_theta_K; }
    const VectorXd& get_theta_K() const {return theta_K; }

    SparseMatrix<double, 0, int>& getA()  {return A;}

    const VectorXd& getW() const {
        return W;
    }

    const VectorXd& getPrevW()  const {
        return prevW;
    }

    void setW(const VectorXd& newW) {
        if (!fix_flag[latent_fix_W]) {
            prevW = W; W = newW;
        }
    }

    void setPrevW(const VectorXd& W) {
        prevW = W;
    }

    // used in block model
    VectorXd getMean() const {
        return mu.cwiseProduct(V - h);
    }

    /*  2 Variance component   */
    const VectorXd& getV() const {
        return V;
    }
    const VectorXd& getPrevV() const {
        return prevV;
    }
    const VectorXd getSV() const {
        return sigma.array().pow(2).matrix().cwiseProduct(V);
    }
    void setPrevV(const VectorXd& V) {
        prevV = V;
    }

    void sample_V(bool posterior=true);

    /*  3 Operator component   */
    const SparseMatrix<double, 0, int>& getK()  {
        return ope->getK();
    }

    /* 4 for optimizer */
    const VectorXd get_parameter() const;
    const VectorXd get_grad();
    void           set_parameter(const VectorXd&);
    void           finishOpt(int i) {fix_flag[i] = 0; }

    // used for general case
    // double function_K(VectorXd& parameter);
    double function_K(const SparseMatrix<double>& K);

    VectorXd grad_theta_K();
    VectorXd grad_theta_mu();
    VectorXd grad_theta_sigma();
    VectorXd grad_theta_sigma_normal(); // grad of sig. only for normal noise

    VectorXd grad_theta_nu() {
        VectorXd grad(n_nu);
        if (n_nu == 1)
            grad(0) = NoiseUtil::grad_theta_nu(noise_type[0], nu[0], V, prevV, h);
        else {
            // for bivaraite case
            int n = V_size / 2;
            grad(0) = NoiseUtil::grad_theta_nu(noise_type[0], nu[0], V.segment(0, n), prevV.segment(0, n), h.segment(0, n));
            grad(1) = NoiseUtil::grad_theta_nu(noise_type[1], nu[1], V.segment(n, n), prevV.segment(n, n), h.segment(n, n));
        }
        return grad;
    }

    Rcpp::List output() const;
};

#endif