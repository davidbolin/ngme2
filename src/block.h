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
#include "var.h"
#include "latent.h"

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using std::vector;

// ---- Structure for random effects ----
// U|V ~ N(0, Sigma)
class Randeff {
  private:
    std::mt19937 randeff_rng;
    string effect_type;

    MatrixXd B_reff, Sigma, invSigma;
    int n_reff, n_cov_params, n_params; // n_params = n_cov_params + 1
    VectorXd U, mu;
    MatrixXd Dd; // duplicated Matrix Dd * vech = vec
    Var var;
  public:
    Randeff(const Rcpp::List& R_randeff, unsigned long seed);
    int get_n_params() const {return n_params;}
    int get_n_reff() const {return n_reff;}
    MatrixXd& get_B_reff() {return B_reff;}
    const VectorXd& get_U() const {return U;}
    void set_U(const VectorXd& U_) {U = U_;}
    VectorXd getMean() const {double V = var.getV()(0); return (V-1)*mu;}

    VectorXd get_parameter();
    VectorXd precond_grad();
    void set_parameter(const VectorXd& p);
};


const int BLOCK_FIX_FLAG_SIZE = 3;
enum Block_fix_flag {
    block_fix_beta, block_fix_theta_mu, block_fix_theta_sigma
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

    // Random effects
    vector<std::unique_ptr<Randeff>> randeffs;

    MatrixXd B_sigma;
    VectorXd noise_sigma, theta_sigma;
    int n_theta_sigma;

    int n_latent; // how mnay latent model
    int n_obs;  // how many observation
    int n_params, n_la_params, n_re_params, n_feff, n_merr, n_reffs;  // number of total params, la params, ...

    // fix estimation
    bool fix_flag[BLOCK_FIX_FLAG_SIZE] {0};

    // controls
    int n_gibbs;
    bool debug, reduce_var;
    double reduce_power, threshold;

    SparseMatrix<double> A, K, G;      // not used: dK, d2K; G = [B_reff A]
    MatrixXd B_reffs;

    vector<std::unique_ptr<Latent>> latents;
    Var var;

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

    int get_n_obs() const {return n_obs;}
    void sampleW_VY();
    void sampleV_WY() {
      if(n_latent > 0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_cond_V();
        }
      }
    }
    void sample_V() {
      if(n_latent > 0){
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_V();
        }
      }
    }
    void setW(const VectorXd&);
    void setPrevW(const VectorXd&);
    void setPrevV(const VectorXd&);

    /* Optimizer related */
    int                  get_n_params() const {return n_params;}

    VectorXd             get_parameter_no_reff() const;
    VectorXd             get_parameter_reff() const;
    void                 set_parameter_no_reff(const VectorXd&);
    void                 set_parameter_reff(const VectorXd&);
    VectorXd             precond_grad_no_reff();
    VectorXd             precond_grad_reff();

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
    VectorXd getMean_reff() const {
        VectorXd mean (n_reffs);
        int pos = 0;
        for (vector<std::unique_ptr<Randeff>>::const_iterator it = randeffs.begin(); it != randeffs.end(); it++) {
            int size = (*it)->get_n_reff();
            mean.segment(pos, size) = (*it)->getMean();
            pos += size;
        }
        return mean;
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
      if (n_latent > 0){
        return Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + var.getV()).cwiseProduct(noise_mu);
      } else {
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

    // for updating hessian
    vector<VectorXd> get_VW() const {
        vector<VectorXd> ret (3);
        ret[0] = var.getV();
        ret[1] = getV();
        ret[2] = getW();
        return ret;
    }

    void set_prev_VW(const vector<VectorXd>& VW) {
        var.setPrevV(VW[0]);
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
