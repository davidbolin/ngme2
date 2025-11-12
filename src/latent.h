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
#include <limits>
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
    latent_fix_theta_nu,
    latent_fix_theta_sigma,
    latent_fix_theta_sigma_normal
};
const int LATENT_FIX_FLAG_SIZE = 7;

class Latent {
protected:
    std::mt19937 latent_rng;
    string model_type;
    vector<string> noise_type;
    bool debug;
    int n_noise, W_size, V_size, n_params, n_var {1}; // n_params=n_theta_K + n_theta_mu + n_theta_sigma + n_var

    // operator (for compute K, for compute numerical gradient, for preconditioner)
    std::shared_ptr<Operator> ope, ope_precond, ope_addeps;
    VectorXd h, theta_K;
    int n_theta_K;
    bool symmetricK, zero_trace;

    struct DerivativeCache {
        bool initialized {false};
        bool rao_blackwell {false};
        VectorXd grad_theta_K;
        VectorXd grad_theta_mu;
        VectorXd grad_theta_sigma;
        VectorXd grad_theta_nu;
        bool grad_K_ready {false};
        bool grad_mu_ready {false};
        bool grad_sigma_ready {false};
        bool grad_nu_ready {false};
    } deriv_cache;

    struct HessianCache {
        bool ready {false};
        // Diagonal blocks
        MatrixXd H_K;      // n_theta_K x n_theta_K
        MatrixXd H_mu;     // n_theta_mu x n_theta_mu
        MatrixXd H_sigma;  // n_theta_sigma x n_theta_sigma
        MatrixXd H_nu;     // n_theta_nu x n_theta_nu (no cross terms with others)
        // Cross blocks (upper-right by convention)
        MatrixXd H_K_mu;       // n_theta_K x n_theta_mu
        MatrixXd H_K_sigma;    // n_theta_K x n_theta_sigma
        MatrixXd H_mu_sigma;   // n_theta_mu x n_theta_sigma
    } hess_cache;

    vector<double> trace;
    VectorXd rb_trace_K, rb_trace_sigma;
    double eps {1e-5};

    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0}, numer_grad {false}, use_iterative_solver {false}, use_same_V {false};
    
    vector<bool> fix_theta_sigma_vec;  // Vector-based fixing for theta_sigma parameters

    vector<std::shared_ptr<Operator>> ope_add_eps;
    vector<SparseMatrix<double>> num_dK;
    // Observation score s = A^T D r passed from Block for Z chain rule
    VectorXd obs_score;

    // mu and sigma, and sigma_normal (special case when using nig_normal case)
    MatrixXd B_mu, B_sigma, B_nu;
    VectorXd theta_mu, theta_sigma, theta_nu;

    // mu = Bmu * theta_mu
    // sigma = exp(Bsigma * theta_sigma)
    // nu = exp(Bnu * theta_nu)
    VectorXd mu, sigma, nu;
    int n_theta_mu, n_theta_sigma, n_theta_nu;

    double nu_lower_bound {1e-3};

    // for numerical gradient and observation mapping
    VectorXd W, prevW, cond_W, V, prevV;
    SparseMatrix<double,0,int> A;
    

    int dim {1}; // noise dimension
    VectorXd p_vec, a_vec, b_vec;
    VectorXd p_inc, a_inc, b_inc;
    bool single_V {false}, share_V {false};

    // Solver controls propagated from control_opt via Block/Latent constructor
    int solver_type_ {0};
    int n_trace_iter_ {8};

    // Short-lived cache for K * W within a Gibbs pass
    Eigen::VectorXd KW_cache_;
    bool KW_valid_ {false};

    // priors
    string prior_K_type, prior_mu_type, prior_sigma_type, prior_nu_type;
    VectorXd prior_K_param, prior_mu_param, prior_sigma_param, prior_nu_param;

    int iter_solver_iter {0};

public:
    Latent(const Rcpp::List&, unsigned long seed);
    ~Latent() {}

    /*  1 Model itself   */
    int get_W_size() const             {return W_size; }
    int get_V_size() const             {return V_size; }
    int get_n_params() const           {return n_params; }
    int get_n_theta_K() const          {return n_theta_K; }
    int get_n_theta_sigma() const      {return n_theta_sigma; }
    
    vector<bool> get_theta_unfixed_sigma() const {return fix_theta_sigma_vec; }

    const VectorXd& get_theta_K() const {return theta_K; }

    SparseMatrix<double, 0, int>& getA()  {return A;}
    const SparseMatrix<double, 0, int>& getZ() const { return ope->getZ(); }
    // Set observation score for Z-chain gradient
    void set_obs_score(const VectorXd& s) { obs_score = s; }

    const VectorXd& getW() const {
        return W;
    }

    const VectorXd& getPrevW()  const {
        return prevW;
    }

    void setW(const VectorXd& newW) {
        if (!fix_flag[latent_fix_W]) { W = newW; }
        if (use_same_V) { prevW = W; } // always update prevW
        invalidate_derivatives();
        invalidate_KW_cache();
    }

    const VectorXd& get_cond_W() const { return cond_W; }
    void set_cond_W(const VectorXd& W) { cond_W = W; invalidate_derivatives(); }
    void setPrevW(const VectorXd& W) { prevW = W; invalidate_derivatives(); }

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
        invalidate_derivatives();
    }
    
    void update_each_iter();
    void sample_cond_V();
    void sample_uncond_V();

    void invalidate_derivatives();
    void update_derivatives(
        bool need_grad_theta_K,
        bool need_grad_theta_mu,
        bool need_grad_theta_sigma,
        bool need_grad_theta_nu,
        bool rao_blackwell
    );

    void compute_theta_K(bool need_grad, bool rao_blackwell);
    void compute_theta_mu(bool need_grad, bool rao_blackwell);
    void compute_theta_sigma(bool need_grad, bool rao_blackwell);
    void compute_theta_nu(bool need_grad);

    // Compute analytic Hessian blocks (skeleton, placeholder implementation for now)
    void compute_hessian_blocks(bool rao_blackwell);

    // Preconditioner using stored Gibbs samples
    MatrixXd preconditioner();

    // Cached K*W accessor (invalidated when W or K changes)
    const Eigen::VectorXd& get_KW_cached();
    inline void invalidate_KW_cache() { KW_valid_ = false; }

    /*  3 Operator component   */
    const SparseMatrix<double, 0, int>& getK()  {
        return ope->getK();
    }

    const SparseMatrix<double, 0, int>& get_dK(int i) {
        return ope->get_dK(i);
    }

    const VectorXd get_BSigma_col(int i)  {
        return B_sigma.col(i);
    }

    // Expose dZ and compute numerical dZ via operator
    const SparseMatrix<double,0,int>& get_dZ(int i) { return ope->get_dZ(i); }

    void set_rb_trace(const VectorXd& rb_trace_K, const VectorXd& rb_trace_sigma) {
        this->rb_trace_K = rb_trace_K;
        this->rb_trace_sigma = rb_trace_sigma;
    }

    /* 4 for optimizer */
    const VectorXd get_parameter();
    const VectorXd get_grad(bool rao_blackwell=FALSE);
    void           set_parameter(const VectorXd&, bool update_dK=FALSE);
    void           finishOpt(int i) {fix_flag[i] = 0; }


    Rcpp::List output() const;
};

#endif
