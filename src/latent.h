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
#include "sample_rGIG.h"

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
const int LATENT_FIX_FLAG_SIZE = 6;

class Latent {
protected:
    std::mt19937 latent_rng;
    string model_type;
    vector<string> noise_type;
    bool debug;
    int n_noise, W_size, V_size, n_params, n_var {1}; // n_params=n_theta_K + n_theta_mu + n_theta_sigma + n_var
    Eigen::VectorXi fix_parameters;

    // operator (for compute K, for compute numerical gradient, for preconditioner)
    std::shared_ptr<Operator> ope, ope_precond, ope_addeps;
    VectorXd h, theta_K;
    int n_theta_K;
    bool symmetricK, zero_trace, use_num_dK {false};

    vector<double> trace;
    VectorXd rb_trace_K, rb_trace_sigma;
    double eps {1e-5};

    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0};
    bool numer_grad {false};

    vector<std::shared_ptr<Operator>> ope_add_eps;
    vector<SparseMatrix<double>> num_dK;

    // mu and sigma, and sigma_normal (special case when using nig_normal case)
    MatrixXd B_mu, B_sigma, B_sigma_normal, B_nu;
    VectorXd theta_mu, theta_sigma, theta_sigma_normal, theta_nu;

    // mu = Bmu * theta_mu
    // sigma = exp(Bsigma * theta_sigma)
    // nu = exp(Bnu * theta_nu)
    VectorXd mu, sigma, sigma_normal, nu;
    int n_theta_mu, n_theta_sigma, n_theta_nu, n_theta_sigma_normal;

    // for numerical gradient.
    VectorXd W, prevW, cond_W, V, prevV;
    SparseMatrix<double,0,int> A;

    int dim {1}; // noise dimension
    VectorXd p_vec, a_vec, b_vec;
    VectorXd p_inc, a_inc, b_inc;
    bool single_V {false}, share_V {false};

    // solver
    cholesky_solver chol_solver_K;
    lu_sparse_solver lu_solver_K;
    cholesky_solver solver_Q; // Q = KT diag(1/SV) K

    // priors
    string prior_K_type, prior_mu_type, prior_sigma_type, prior_nu_type;
    VectorXd prior_K_param, prior_mu_param, prior_sigma_param, prior_nu_param;
public:
    Latent(const Rcpp::List&, unsigned long seed);
    ~Latent() {}

    /*  1 Model itself   */
    int get_W_size() const             {return W_size; }
    int get_V_size() const             {return V_size; }
    int get_n_params() const           {return n_params; }
    int get_n_theta_K() const          {return n_theta_K; }
    int get_n_theta_sigma() const      {return n_theta_sigma; }

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

    const VectorXd& get_cond_W() const { return cond_W; }
    void set_cond_W(const VectorXd& W) { cond_W = W; }
    void setPrevW(const VectorXd& W) { prevW = W; }

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
    const VectorXd getPrevSV() const {
        return sigma.array().pow(2).matrix().cwiseProduct(prevV);
    }
    void setPrevV(const VectorXd& V) {
        prevV = V;
    }

    void update_each_iter(bool initialization=false, bool update_dK=false);
    void sample_cond_V();
    void sample_uncond_V();

    // pi(W|V)
    double logd_W_given_V(const VectorXd& W, const SparseMatrix<double>& K, const VectorXd& mu, const VectorXd& sigma, const VectorXd& V);
    // pi(KW|V) fix K
    double logd_KW_given_V(const VectorXd& mu, const VectorXd& sigma, const VectorXd& V);

    // pi(W|V) * pi(V)
    double log_density(const VectorXd& parameter, bool precond_K = false);
    MatrixXd precond(bool precond_K = false, double eps = 1e-5);

    /*  3 Operator component   */
    const SparseMatrix<double, 0, int>& getK()  {
        return ope->getK();
    }

    const SparseMatrix<double, 0, int>& get_dK(int i);

    const VectorXd get_BSigma_col(int i)  {
        return B_sigma.col(i);
    }

    void set_rb_trace(const VectorXd& rb_trace_K, const VectorXd& rb_trace_sigma) {
        this->rb_trace_K = rb_trace_K;
        this->rb_trace_sigma = rb_trace_sigma;
    }

    /* 4 for optimizer */
    const VectorXd get_parameter();
    const VectorXd get_grad(bool rao_blackwell=FALSE);
    void           set_parameter(const VectorXd&, bool update_dK=FALSE);
    void           finishOpt(int i) {fix_flag[i] = 0; }

    VectorXd grad_theta_K(bool rao_blackwell=FALSE);
    VectorXd grad_theta_mu(bool rao_blackwell=FALSE);
    VectorXd grad_theta_sigma(bool rao_blackwell=FALSE);
    VectorXd grad_theta_sigma_normal(bool rao_blackwell=FALSE); // grad of sig. only for normal noise
    VectorXd grad_theta_nu();

    Rcpp::List output() const;
};

#endif