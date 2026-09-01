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
    std::shared_ptr<Operator> ope;
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
        MatrixXd H_nu;     // n_theta_nu x n_theta_nu
        // Cross blocks (upper-right by convention)
        MatrixXd H_K_mu;       // n_theta_K x n_theta_mu
        MatrixXd H_K_sigma;    // n_theta_K x n_theta_sigma
        MatrixXd H_mu_sigma;   // n_theta_mu x n_theta_sigma
    } hess_cache;

    // Cached gradient and preconditioner data
    VectorXd last_gradient_;
    bool     grad_cache_valid_ {false};
    bool     grad_cache_rb_mode_ {false};
    MatrixXd last_precond_;
    bool     precond_cache_valid_ {false};
    bool     state_ready_ {false};
    bool     state_has_precond_terms_ {false};

    vector<double> trace;
    double eps {1e-5};

    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0}, numer_grad {false}, use_iterative_solver {false}, use_same_V {false};

    vector<bool> fix_theta_sigma_vec;  // Vector-based fixing for theta_sigma parameters

    // mu and sigma, and sigma_normal (special case when using nig_normal case)
    MatrixXd B_mu, B_sigma, B_nu;
    VectorXd theta_mu, theta_sigma, theta_nu;

    // mu = Bmu * theta_mu
    // sigma = exp(Bsigma * theta_sigma)
    // nu = nu_lower_bound + exp(Bnu * theta_nu)
    VectorXd mu, sigma, nu;
    int n_theta_mu, n_theta_sigma, n_theta_nu;

  double nu_lower_bound {0.0};

    // for numerical gradient and observation mapping
    VectorXd W, prevW, cond_W, V, prevV;
    SparseMatrix<double,0,int> A;

    int dim {1}; // noise dimension
    VectorXd p_vec, a_vec, b_vec;
    VectorXd p_inc, a_inc, b_inc;
    bool single_V {false}, share_V {false};

    // Solver controls propagated from control_opt via Block/Latent constructor
    int solver_type_ {0};
    int nonsym_solver_ {0};
    int n_trace_iter_ {8};
    bool robust_ {false};
    // Standardised (Cabral, Bolin & Rue 2023, sec 2.1) coordinates for the
    // stationary NIG noise. The optimiser works in
    //     t = (log sigma_marg, zeta, log eta),
    //     eta = 1/nu,  zeta = mu/sigma,  sigma_marg = sqrt(sigma^2 + mu^2 eta),
    // where sigma_marg is the marginal SD, since Var = h (sigma^2 + mu^2/nu).
    // In the native coordinates that variance is shared by all three
    // parameters, so the likelihood has a long flat ridge along
    // sigma^2 + mu^2/nu = const; here that ridge is the sigma_marg axis.
    // 0 = native, 1 = standardised, 2 = additionally
    // orthogonalised zeta* = zeta sqrt(eta), eta* = eta / xi^2 with
    // xi = 1 + zeta*^2 - |zeta*| sqrt(1 + zeta*^2), which makes the large-
    // deviation rate (hence the kurtosis) invariant to skewness.
    int nig_param_mode_ {0};
    bool nig_std_active() const;
    VectorXd nig_std_from_native() const;
    // t -> (theta_mu, theta_sigma, theta_nu)
    VectorXd nig_native_from_t(const VectorXd &t) const;
    void nig_std_to_native(const VectorXd &t, double &theta_mu_out,
                           double &theta_sigma_out, double &theta_nu_out) const;
    // d(native)/d(t), by central differences on nig_native_from_t. The map is a
    // handful of flops, so this is far cheaper than the likelihood and avoids a
    // second hand-derivation for each mode.
    MatrixXd nig_std_jacobian() const;

    // priors
    std::vector<string> prior_K_type;
    std::vector<VectorXd> prior_K_param;
    std::vector<string> prior_K_target;
    string prior_mu_type, prior_sigma_type, prior_nu_type;
    VectorXd prior_mu_param, prior_sigma_param, prior_nu_param;
    string prior_mu_target {"coef"};
    string prior_sigma_target {"coef"};
    string prior_nu_target {"coef"};

    int iter_solver_iter {0};

public:
    Latent(const Rcpp::List&, unsigned long seed);
    ~Latent() {}

    /*  1 Model itself   */
    int get_W_size() const             {return W_size; }
    int get_V_size() const             {return V_size; }
    int get_n_params() const           {return n_params; }
    // Empty unless the standardised coordinates are active. Block applies this
    // once the whole native gradient (including its own RB and dZ terms) is
    // assembled -- transforming earlier would mix coordinate systems.
    MatrixXd get_nig_std_jacobian() const {
      return nig_std_active() ? nig_std_jacobian() : MatrixXd();
    }
    int get_n_theta_K() const          {return n_theta_K; }
    int get_n_theta_sigma() const      {return n_theta_sigma; }
    int get_n_theta_mu() const         {return n_theta_mu; }
    int get_n_theta_nu() const         {return n_theta_nu; }

    vector<bool> get_theta_unfixed_sigma() const {return fix_theta_sigma_vec; }

    const VectorXd& get_theta_K() const {return theta_K; }

    SparseMatrix<double, 0, int>& getA()  {return A;}
    const SparseMatrix<double, 0, int>& getZ() const { return ope->getZ(); }

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

    void update_each_iter(bool need_precond = false);
    void sample_cond_V();
    void sample_uncond_V();

    // True when a call to sample_cond_V() / sample_uncond_V() can actually
    // change V. For a purely Gaussian (and non-fixed-V) latent both samplers
    // skip every component, so V is left untouched. Used by BlockModel to decide
    // whether QQ has to be reassembled and refactorized between Gibbs draws.
    bool V_may_change() const {
        if (fix_flag[latent_fix_V]) return false;
        // the single_V + share_V branch of sample_cond_V() redraws V without
        // consulting noise_type, so treat it as always changing.
        if (single_V && share_V) return true;
        for (const string& nt : noise_type)
            if (nt != "normal") return true;
        return false;
    }

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
    MatrixXd preconditioner() const;

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
    const SparseMatrix<double,0,int>& get_d2Z(int i, int j) const { return ope->get_d2Z(i,j); }

    /* 4 for optimizer */
    const VectorXd get_parameter();
    const VectorXd get_grad();
    void           compute_grad_and_hessian(bool rao_blackwell, bool with_precond);
    void           set_parameter_and_update(const VectorXd&, bool with_precond);
    void           finishOpt(int i) {fix_flag[i] = 0; }

    Rcpp::List output() const;
};

#endif
