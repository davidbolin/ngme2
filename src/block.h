/*
BlockModel
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
#include "var.h"
#include "latent.h"

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using std::vector;

const int BLOCK_FIX_FLAG_SIZE = 6;

enum Block_fix_flag {
    block_fix_beta, block_fix_theta_mu, block_fix_theta_sigma, block_fix_theta_V,
    block_fix_V
};

class BlockModel : public Model {
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
    int n_obs; // how many observation
    int n_params, n_la_params, n_feff, n_merr;  // number of total params, la params, ...

    // fix estimation
    bool fix_flag[BLOCK_FIX_FLAG_SIZE] {0};

    // controls
    int n_gibbs;
    bool debug,opt_beta, kill_var;
    double kill_power, threshold, termination;

    SparseMatrix<double> A, K;      // not used: dK, d2K;

    // std::unique_ptr<Var> var;
    Var var;

    // debug
    bool fix_merr;
    bool fixblockV;

    // optimize related
    VectorXd stepsizes, gradients;
    int counting {0};
    VectorXd indicate_threshold, steps_to_threshold;

    // No initializer
    std::vector<std::unique_ptr<Latent>> latents;
    cholesky_solver chol_Q, chol_QQ;
    SparseLU<SparseMatrix<double> > LU_K;

    VectorXd fixedW;

    // record trajectory
    vector<vector<double>> beta_traj;
    vector<vector<double>> theta_mu_traj;
    vector<vector<double>> theta_sigma_traj;
    vector<double>   theta_V_traj;
public:
    // BlockModel() {}
    BlockModel(Rcpp::List& block_model, unsigned long seed);
    virtual ~BlockModel() {}

    /* Gibbs Sampler */
    void burn_in(int iterations) {
        for (int i=0; i < iterations; i++) {
            sampleW_VY();
            sampleV_WY();
            sample_cond_block_V();
        }
      if (debug) std::cout << "Finish burn in period." << std::endl;
    }

    void sampleW_VY();
    void sampleV_WY() {
      if(n_latent >0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_cond_V();
        }
      }
    }
    void sample_V() {
      if(n_latent >0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_V();
        }
      }
    }
    void setW(const VectorXd&);

    /* Optimizer related */
    VectorXd             get_parameter() const;
    VectorXd             get_stepsizes() const {return stepsizes;}
    void                 set_parameter(const VectorXd&);
    VectorXd             grad();
    SparseMatrix<double> precond() const;

    void                 examine_gradient();
    void                 sampleW_V();

    /* Aseemble */
    void assemble() {
        int nrow = 0;
        int ncol = 0;
        for (std::vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
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
        for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            mean.segment(pos, size) = (*it)->getMean();
            pos += size;
        }
        return mean;
    }

    VectorXd getV() const {
        VectorXd V (V_sizes);
        int pos = 0;
        for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            V.segment(pos, size) = (*it)->getV();
            pos += size;
        }

        return V;
    }

    // return sigma * V
    VectorXd getSV() const {
        VectorXd SV (V_sizes);
        int pos = 0;
        for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_V_size();
            SV.segment(pos, size) = (*it)->getSV();
            pos += size;
        }

        return SV;
    }

    VectorXd getW() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->getW();
            pos += size;
        }
        return W;
    }

    VectorXd getPrevW() const {
        VectorXd W (W_sizes);
        int pos = 0;
        for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->get_W_size();
            W.segment(pos, size) = (*it)->getPrevW();
            pos += size;
        }
        return W;
    }

    VectorXd get_residual() const {
      if(n_latent>0){
        return Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + var.getV()).cwiseProduct(noise_mu);
      }else{
        return Y  - X * beta - (-VectorXd::Ones(n_obs) + var.getV()).cwiseProduct(noise_mu);
      }
    }

    void sample_cond_block_V() {
        if (family == "nig") {
            VectorXd residual = get_residual();
            VectorXd a_inc_vec = noise_mu.cwiseQuotient(noise_sigma).array().pow(2);
            VectorXd b_inc_vec = (residual + var.getV().cwiseProduct(noise_mu)).cwiseQuotient(noise_sigma).array().pow(2);
            var.sample_cond_V(a_inc_vec, b_inc_vec);
        }
    }

    // --------- Fixed effects and Measurement error  ------------
    VectorXd grad_beta();

    VectorXd get_theta_merr() const;
    VectorXd grad_theta_mu();
    VectorXd grad_theta_sigma();
    VectorXd grad_theta_merr();
    void set_theta_merr(const VectorXd& theta_merr);

    // get length of W,V of iterations
    Rcpp::List sampling(int iterations, bool posterior);

    // return output
    Rcpp::List output() const;
};

// ---- inherited functions ------
/* the way structuring the parameter
    latents[1].get_parameter
        ...
    beta (fixed effects)
    measurement noise
*/

// Not Implemented
inline SparseMatrix<double> BlockModel::precond() const {
    SparseMatrix<double> precond;
    std::cout << "Not implemented \n";
    throw;
    return precond;
}

#endif
