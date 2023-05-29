/*
BlockModel - for ngme each replicate
*/
#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

// #define EIGEN_USE_MKL_ALL
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>

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

const int BLOCK_FIX_FLAG_SIZE = 5;
enum Block_fix_flag {
    block_fix_beta, block_fix_theta_mu, block_fix_theta_sigma, blcok_fix_V, block_fix_nu
};

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
    VectorXd beta;
    MatrixXd B_mu;
    VectorXd noise_mu, theta_mu;
    int n_theta_mu;

    MatrixXd B_sigma;
    VectorXd noise_sigma, theta_sigma;
    int n_theta_sigma;

    int n_latent; // how mnay latent model
    int n_obs;  // how many observation
    int n_la_params, n_feff, n_merr, n_repl; // number of total replicates(blocks)

    // Correlated measurement error
    bool corr_measure;
    double rho {0};
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

    SparseMatrix<double> A, K, Q, QQ;      // not used: dK, d2K; G = [B_reff A]
    SimplicialLLT<SparseMatrix<double, Lower>> Q_eps_solver;

    vector<std::unique_ptr<Latent>> latents;
    VectorXd p_vec, a_vec, b_vec, noise_V, noise_prevV;
    double nu {1};

    // optimize related
    VectorXd stepsizes, gradients;
    int counting {0};
    VectorXd indicate_threshold, steps_to_threshold;
    int curr_iter; // how many times set is called.

    // solvers
    cholesky_solver chol_Q, chol_QQ;
    SparseLU<SparseMatrix<double>> LU_K;

    std::string par_string;
public:
    // BlockModel() {}
    BlockModel(const Rcpp::List& block_model, unsigned long seed);
    ~BlockModel() = default;

    /* Gibbs Sampler */
    void burn_in(int iterations) {
        for (int i=0; i < iterations; i++) {
            sampleW_VY();
            sampleV_WY();
            sample_noise_V();
        }
    }

    int get_n_obs() const {return n_obs;}
    void sampleW_VY();
    void sampleV_WY() {
      if(n_latent > 0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_V(true);
        }
      }
    }

    void sample_uncond_V() {
      if(n_latent > 0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_V(false);
        }
      }
    }

    void setW(const VectorXd&);
    void setPrevW(const VectorXd&);
    void setPrevV(const VectorXd&);

    /* Optimizer related */
    int                  get_n_params() const {return n_params;}

    VectorXd             get_parameter() const;
    void                 set_parameter(const VectorXd&);
    VectorXd             precond_grad();

    MatrixXd             precond() const;
    // VectorXd             grad() {return precond_grad();}

    void                 examine_gradient();
    void                 sampleW_V();

    /* Aseemble */
    void assemble() {
        int nrow = 0;
        int ncol = 0;
        for (vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
            setSparseBlock(&K,   nrow, ncol, (*it)->getK());
            // setSparseBlock(&dK,  n, n, (*it)->get_dK());
            // setSparseBlock(&d2K, n, n, (*it)->get_d2K());
            nrow += (*it)->get_V_size();
            ncol += (*it)->get_W_size();
        }
    }

    // return mean = mu*(V-h)
    VectorXd getMean() const {
        VectorXd mean (V_sizes);
        int pos = 0;
        for (vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            mean.segment(pos, size) = (*it)->getMean();
            pos += size;
        }
        return mean;
    }

    VectorXd getV() const {
        VectorXd V (V_sizes);
        int pos = 0;
        for (vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            V.segment(pos, size) = (*it)->getV();
            pos += size;
        }

        return V;
    }

    VectorXd getPrevV() const {
        VectorXd V (V_sizes);
        int pos = 0;
        for (vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
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
        for (vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            SV.segment(pos, size) = (*it)->getSV();
            pos += size;
        }

        return SV;
    }

    VectorXd getW() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->getW();
            pos += size;
        }
        return W;
    }

    VectorXd getPrevW() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->getPrevW();
            pos += size;
        }
        return W;
    }

    VectorXd get_residual() const {
      if (n_latent > 0) {
        return Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
      } else {
        return Y  - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
      }
    }

    // residual_part = residual + AW
    VectorXd get_residual_part() const {
        return Y - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
    }

    void sample_noise_V(bool posterior = true) {
        if (family == "normal" || fix_flag[blcok_fix_V]) return;
        noise_prevV = noise_V;

        if (posterior) {
            VectorXd a_inc_vec = noise_mu.cwiseQuotient(noise_sigma).array().pow(2);
            VectorXd b_inc_vec = (get_residual() + noise_V.cwiseProduct(noise_mu)).cwiseQuotient(noise_sigma).array().pow(2);
            double dim = 1;
            VectorXd p_vec_new = p_vec - VectorXd::Constant(n_obs, 0.5 * dim);
            VectorXd a_vec_new = a_vec + a_inc_vec;
            VectorXd b_vec_new = b_vec + b_inc_vec;
            // noise_V = rGIG_cpp(p_vec_new, a_vec_new, b_vec_new, rng());
            NoiseUtil::sample_V(noise_V, family, p_vec_new, a_vec_new, b_vec_new, rng);
        } else {
            // noise_V = rGIG_cpp(p_vec, a_vec, b_vec, rng());
            NoiseUtil::sample_V(noise_V, family, p_vec, a_vec, b_vec, rng);
        }
    }

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
    Rcpp::List sampling(int n, bool posterior);
    Rcpp::List output() const;
    std::string get_par_string() const {return par_string;}

    double th2rho(double th) const {return (-1 + 2*exp(th) / (1+exp(th)));}
    double rho2th(double r) const {return (log((-1-r)/(-1+r)));}

    // drho / dtheta
    // double dtheta_rho(double th) const {return 2 * exp(th) / pow(1+exp(th), 2);}
    double dtheta_th(double rho) const {return (1-rho*rho) / 2;}
};

// ---- inherited functions ------
/* the way structuring the parameter
    latents[1].get_parameter
        ...
    beta (fixed effects)
    measurement noise
*/

// Not Implemented
inline MatrixXd BlockModel::precond() const {
    MatrixXd precond;
    std::cout << "Not implemented \n";
    throw;
    return precond;
}


#endif
