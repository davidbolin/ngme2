/*
BlockModel
*/

#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>

#include <string>
#include <vector>
#include <iostream>

#include "include/timer.h"
#include "include/solver.h"
#include "include/MatrixAlgebra.h"
#include "model.h"
#include "latent.h"

#include "latents/ar1.h"

using Eigen::SparseMatrix;

const int latent_para = 4;

class BlockModel : public Model {
protected:

// n_regs   = row(A1) + ... + row(An)
// n_paras = 4 * n_latent + 1
    int n_obs, n_latent, n_regs, n_paras, n_gibbs;
    std::vector<Latent*> latents;

    MatrixXd X;
    VectorXd Y, beta;
    
    string family;
    bool opt_fix_effect {true};
    double sigma_eps;
    
    SparseMatrix<double> A, K, dK, d2K;
    cholesky_solver chol_Q;
    SparseLU<SparseMatrix<double> > LU_K;
public:
    BlockModel() {}

    BlockModel(MatrixXd X,
               VectorXd Y, 
               string family,
               int n_regs,
               Rcpp::List latents_in,
               int n_gibbs,
               VectorXd beta,
               double sigma_eps,
               bool opt_fix_effect) 
    : n_obs(Y.size()),
      n_latent(latents_in.size()), 
      n_regs(n_regs),
      n_paras(n_latent * latent_para + 1),
      n_gibbs(n_gibbs),
      
      X(X),
      Y(Y), 
      beta(beta),

      family(family),
      opt_fix_effect(opt_fix_effect),
      sigma_eps(sigma_eps), 

      A(n_obs, n_regs), 
      K(n_regs, n_regs), 
      dK(n_regs, n_regs),
      d2K(n_regs, n_regs)
    {
    // Init each latent model
    for (int i=0; i < n_latent; ++i) {
        Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);

        // construct acoording to models
        string type = latent_in["type"];
        if (type == "ar1") {
            latents.push_back(new AR(latent_in) );
        }
    }
    
    /* Fixed effects */
    if (opt_fix_effect) {
        int n_beta = beta.size();
        n_paras = latent_para + 1 + n_beta;
    }

    /* Init variables: h, A */
    int n = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        setSparseBlock(&A,   0, n, (*it)->getA());            
        n += (*it)->getSize();
    }
    assemble();

        VectorXd inv_SV = VectorXd::Constant(n_regs, 1).cwiseQuotient(getSV());
        SparseMatrix<double> QQ = K.transpose() * inv_SV.asDiagonal() * K + pow(sigma_eps, 2) * A.transpose() * A;
    
    chol_Q.analyze(QQ);
    LU_K.analyzePattern(K);

    // sample W
    
    sampleW_VY();
    // sampleW_V();
    }

    void sampleW_VY();
    void sampleW_V();

    void sampleV_WY() {
        for (unsigned i=0; i < n_latent; i++) {
            (*latents[i]).sample_cond_V();
        }
    }

    void setW(const VectorXd&);
    VectorXd             get_parameter() const;
    void                 set_parameter(const VectorXd&);
    VectorXd             grad();
    SparseMatrix<double> precond() const;

    void assemble() {
        int n = 0;
        for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
            setSparseBlock(&K,   n, n, (*it)->getK());      
            setSparseBlock(&dK,  n, n, (*it)->get_dK());   
            setSparseBlock(&d2K, n, n, (*it)->get_d2K()); 
            
            n += (*it)->getSize();
        }
    }

    // return mean
    VectorXd getMean() const {
        VectorXd mean (n_regs);
        int pos = 0;
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->getSize();
            mean.segment(pos, size) = (*it)->getMean();
            pos += size;
        }
        return mean;
    }

    // return sigma * V
    VectorXd getSV() const {
        VectorXd SV (n_regs);
        int pos = 0;
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->getSize();
            SV.segment(pos, size) = (*it)->getSV();
            pos += size;
        }

        return SV;
    }

    VectorXd getW() const {
        VectorXd W (n_regs);
        int pos = 0;
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->getSize();
            W.segment(pos, size) = (*it)->getW();
            pos += size;
        }
        return W;
    }

    double get_theta_sigma_eps() const {
        return log(sigma_eps);
        // return sigma_eps;
    }

    double grad_theta_sigma_eps() const {
        double g = 0;
        if (family=="normal") {
            VectorXd W = getW();
            VectorXd tmp = Y - A * W - X * beta;
            double norm2 =  tmp.dot(tmp);
            
            g = -(1.0 / sigma_eps) * n_obs  + pow(sigma_eps, -3) * norm2;
            
            double hess = -2.0 * n_obs * pow(sigma_eps, -2);
            // double hess = 1.0 * n_obs * pow(sigma_eps, -2)  - 3 * pow(sigma_eps, -4) * norm2;

            g = g / (hess * sigma_eps + 2 * g); 
        }
        return g;
    }

    void set_theta_sgima_eps(double theta) {
        sigma_eps = exp(theta);
        // sigma_eps = theta;
    }

// Q: it's converging too slow
// fixed effects
    VectorXd grad_beta() const {
        MatrixXd hess = X.transpose() * X;
        VectorXd grads = X.transpose() * (Y - X*beta - A*getW());
std::cout << "old grads=" << -grads << std::endl;
        grads = hess.ldlt().solve(grads);
std::cout << "new grads=" << -grads << std::endl;        
        return -grads;
    }

// FOR TESTING
Rcpp::List testResult() {
Rcpp::List res;
A.makeCompressed();
K.makeCompressed();
dK.makeCompressed();
d2K.makeCompressed();
res["A"] = A;
res["K"] = K;
res["dK"] = dK;
res["d2K"] = d2K;
res["V"] = getSV();
res["W"] = getW();
return res;
}

};

// ---- inherited functions ------

inline VectorXd BlockModel::get_parameter() const {
    VectorXd thetas (n_paras);
    int pos = 0;
    for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
        VectorXd theta = (*it)->getTheta();
        thetas.segment(pos, theta.size()) = theta;
        pos += theta.size();
    }
    
    // sigma_eps
    thetas(n_paras-1) = get_theta_sigma_eps();
    
    // fixed effects
    if (opt_fix_effect) {
        int n_beta = beta.size();
        thetas.segment(n_paras - n_beta-1, n_beta) = beta;
    }

    return thetas;
}

inline VectorXd BlockModel::grad() {
    VectorXd avg_gradient = VectorXd::Zero(n_paras);
    
    for (int i=0; i < n_gibbs; i++) {
        
        // stack grad
        VectorXd gradient = VectorXd::Zero(n_paras);
        
        // get grad for each latent
        int pos = 0;
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int theta_len = (*it)->getThetaSize();
auto timer_computeg = std::chrono::steady_clock::now();
            gradient.segment(pos, theta_len) = (*it)->getGrad();
std::cout << "timer_computeg (ms): " << since(timer_computeg).count() << std::endl;   
            pos += theta_len;
        }

        avg_gradient += gradient;

        // gibbs sampling
        sampleV_WY(); 
auto timer_sampleW = std::chrono::steady_clock::now();
        sampleW_VY();
        // sampleW_V();
std::cout << "sampleW (ms): " << since(timer_sampleW).count() << std::endl;
    }
    avg_gradient = (1.0/n_gibbs) * avg_gradient;

    // sigma_eps 
    avg_gradient(n_paras-1) = grad_theta_sigma_eps();

    // fixed effects
    if (opt_fix_effect) {
        int n_beta = beta.size();
        avg_gradient.segment(n_paras - n_beta-1, n_beta) = grad_beta();
    }

    return avg_gradient;
}

inline void BlockModel::set_parameter(const VectorXd& Theta) {
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        int theta_len = (*it)->getThetaSize();
        VectorXd theta = Theta.segment(pos, theta_len);
        (*it)->setTheta(theta);
        pos += theta_len;
    }
    // sigma_eps
    set_theta_sgima_eps(Theta(n_paras-1));
    
    // fixed effects
    if (opt_fix_effect) {
        int n_beta = beta.size();
        beta = Theta.segment(n_paras - n_beta-1, n_beta);
    }

    assemble(); //update K,dK,d2K after
// std::cout << "Theta=" << Theta <<std::endl;
}


inline SparseMatrix<double> BlockModel::precond() const {
    SparseMatrix<double> precond;

    return precond;
}



#endif