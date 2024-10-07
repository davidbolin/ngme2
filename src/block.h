/*
BlockModel - for ngme each replicate
*/
#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <memory>
#include "include/timer.h"
#include "include/solver.h"
#include "include/MatrixAlgebra.h"
#include "model.h"
#include "latent.h"
#include "noise.h"

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using std::vector;

const int BLOCK_FIX_FLAG_SIZE = 7;
enum Block_fix_flag {
    block_fix_beta, 
    // noise_parameters
    block_fix_theta_mu, block_fix_theta_sigma, block_fix_theta_nu, block_fix_rho, blcok_fix_V, 
    // currently not supported
    block_fix_theta_sigma_normal
};

enum precond_type {precond_none, precond_fast, precond_full};

class BlockModel {
protected:
// W_sizes = row(A1) + ... + row(An)
    // general
    std::mt19937 rng;

    MatrixXd X;
    VectorXd Y;
    int W_sizes, V_sizes; //V_sizes = sum(nrow(K_i))
    string family;

    // Fixed effects and Measurement noise
    VectorXd theta_mu, theta_sigma, theta_nu, rho, beta;

    MatrixXd B_mu, B_sigma, B_nu;
    VectorXd noise_mu, noise_sigma, noise_nu;
    int n_theta_mu, n_theta_sigma, n_theta_nu, n_rho;

    double nu_lower_bound {1e-3};

    int n_latent; // how mnay latent model
    int n_obs;  // how many observation

    // n_merr = noise_nu.size + noise_sigma.size + noise_nu.size + rho.size beta.size
    int n_la_params, n_feff, n_merr, n_repl; // number of total replicates(blocks)

    // Correlated measurement error
    bool corr_measure;
    vector<int> cor_rows, cor_cols;
    vector<bool> has_correlation;
    int n_corr_pairs;
    SparseMatrix<double> Q_eps, dQ_eps;
    int n_params;

    // fix estimation
    bool fix_flag[BLOCK_FIX_FLAG_SIZE] {0};

    // controls
    int n_gibbs;
    bool debug, reduce_var;
    double reduce_power, threshold;

    SparseMatrix<double> A, K, Q, QQ, pmat, pmat_inv;

    vector<std::shared_ptr<Latent>> latents;
    VectorXd p_vec, a_vec, b_vec, noise_V, noise_prevV;
    // double nu {1};

    // optimize related
    VectorXd stepsizes;
    int counting {0};
    VectorXd indicate_threshold, steps_to_threshold;
    int curr_iter; // how many times set is called.

    // solvers
    // SimplicialLLT<SparseMatrix<double, Lower>> Q_eps_solver;
    iterative_solver iterative_QQ;
    cholesky_solver chol_Q, chol_QQ, chol_Q_eps;
    SparseLU<SparseMatrix<double>> LU_K;
    double logdet_Q_eps;

    bool all_gaussian, rao_blackwell, shared_sigma, use_iterative_solver; // No need for gibbs sampling
    std::string par_string;
    VectorXd rb_trace_noise_sigma;

    // priors
    string prior_mu_type, prior_sigma_type, prior_nu_type;
    VectorXd prior_mu_param, prior_sigma_param, prior_nu_param;

    // For computing RB gradient_K
    vector<vector<SparseMatrix<double>>> block_dK;

    // clock
    std::chrono::milliseconds sampling_time {0}, update_time {0};
public:
    // BlockModel() {}
    BlockModel(const Rcpp::List& block_model, unsigned long seed);
    ~BlockModel() = default;

    /* Gibbs Sampler */
    void burn_in(int);

    int get_n_obs() const {return n_obs;}
    void sampleW_VY(bool burn_in = false);
    
    double log_likelihood();

    void sample_cond_V(bool update_Q = true) {
      if(n_latent > 0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_cond_V();
        }
        
        if (update_Q) update_QQ();
      }
    }

    void sample_uncond_V() {
      if(n_latent > 0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_uncond_V();
        }
      }
    }

    void update_QQ();
    void update_Q_eps(double rho);

    void setW(const VectorXd&);
    void setPrevW(const VectorXd&);
    void setPrevV(const VectorXd&);
    void set_cond_W(const VectorXd&);

    /* Optimizer related */
    int                  get_n_params() const {return n_params;}

    VectorXd             get_parameter();
    void                 set_parameter(const VectorXd&);

    VectorXd             grad();
    MatrixXd             precond(int strategy=0, double eps=1e-5);

    void                 examine_gradient();
    void                 sampleW_V();

    /* Aseemble */
    void assemble() {
        int nrow = 0; int ncol = 0;
        for (vector<std::shared_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
            setSparseBlock(&K,   nrow, ncol, (*it)->getK());
            // setSparseBlock(&dK,  n, n, (*it)->get_dK());
            // setSparseBlock(&d2K, n, n, (*it)->get_d2K());
            nrow += (*it)->get_V_size();
            ncol += (*it)->get_W_size();
        }
    }

    // assemble dK, and dK_sigma
    void assemble_dK();

    // tr(QQ^-1 dK^T diag(1/SV) K)
    void compute_rb_trace();

    // return mean = mu*(V-h)
    VectorXd getMean() const {
        VectorXd mean (V_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            mean.segment(pos, size) = (*it)->getMean();
            pos += size;
        }
        return mean;
    }

    VectorXd getV() const {
        VectorXd V (V_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            V.segment(pos, size) = (*it)->getV();
            pos += size;
        }

        return V;
    }

    VectorXd getPrevV() const {
        VectorXd V (V_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            V.segment(pos, size) = (*it)->getPrevV();
            pos += size;
        }

        return V;
    }

    // return sigma * V
    VectorXd getSV() const {
        VectorXd SV (V_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            SV.segment(pos, size) = (*it)->getSV();
            pos += size;
        }

        return SV;
    }

    VectorXd getW() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->getW();
            pos += size;
        }
        return W;
    }

    VectorXd get_cond_W() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->get_cond_W();
            pos += size;
        }
        return W;
    }

    VectorXd getPrevW() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->getPrevW();
            pos += size;
        }
        return W;
    }

    VectorXd get_residual(bool rao_blackwell=false) const {
      if (n_latent > 0 && !rao_blackwell) {
        return Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
      } else if (n_latent > 0 && rao_blackwell) {
        return Y - A * get_cond_W() - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
      } else {
        return Y  - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
      }
    }

    // residual_part = Y - X beta - (1 - V) mu
    VectorXd get_residual_part() const {
        return Y - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
    }

    void sample_cond_noise_V(bool posterior = true);

    // for updating hessian
    vector<VectorXd> get_VW() const {
        vector<VectorXd> ret (3);
        ret[0] = noise_V;
        ret[1] = getV();
        ret[2] = getW();
        return ret;
    }

    void set_prev_VW(const vector<VectorXd>& VW) {
        noise_prevV = VW[0];
        setPrevV(VW[1]);
        setPrevW(VW[2]);
    }


    // --------- Fixed effects and Measurement error  ------------
    VectorXd grad_beta();

    VectorXd get_theta_merr() const;
    VectorXd grad_theta_mu();
    VectorXd grad_theta_sigma();
    VectorXd grad_theta_merr();
    void set_theta_merr(const VectorXd& theta_merr);

    // get length of W,V of iterations
    Rcpp::List sampling(
        int n, 
        int n_burnin,
        bool posterior, 
        const SparseMatrix<double>& A
    );
    
    Rcpp::List sampling(int n, int n_burnin, bool posterior) {
        return sampling(n, n_burnin, posterior, A);
    }

    Rcpp::List output() const;
    std::string get_par_string() const {return par_string;}

    static double th2rho(double th) {return (-1 + 2*exp(th) / (1+exp(th)));}
    static double rho2th(double r) {return (log((-1-r)/(-1+r)));}

    // drho / dtheta
    // double dtheta_rho(double th) const {return 2 * exp(th) / pow(1+exp(th), 2);}
    static double dtheta_th(double rho) {return (1-rho*rho) / 2;}

    // compute numerical hessian for theta = c(beta, theta_mu, theta_sigma, nu)
    MatrixXd num_h_no_latent(const VectorXd& v, double eps = 1e-5);
    double logd_no_latent(const VectorXd& v);

    void test_in_the_end();
};

#endif
