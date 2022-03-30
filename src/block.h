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

/* config */
    unsigned n_gibbs;

// n_regs   = row(A1) + ... + row(An)
// n_paras = 4 * n_latent + 1
    unsigned n_obs, n_latent, n_regs, n_paras;
    std::vector<Latent*> latents;

    string family;
    double sigma_eps;
    
    MatrixXd X;
    VectorXd Y, beta;
    bool opt_fix_effect {true};

    SparseMatrix<double> A, K, dK, d2K;

    cholesky_solver chol_Q;
public:
    BlockModel() {}
    // main constructor
    BlockModel(MatrixXd X,
               VectorXd Y, 
               string family,
               int n_regs,
               Rcpp::List latents_in,
               int n_gibbs,
               VectorXd beta,
               bool opt_fix_effect) 
    : n_gibbs(n_gibbs),
      n_obs(Y.size()),
      n_paras(n_latent * latent_para + 1),
      n_latent(latents_in.size()), 
      n_regs(n_regs),
      sigma_eps(2),
      Y(Y), 
      X(X),
      beta(beta),
      family(family),
      A(n_obs, n_regs), 
      K(n_regs, n_regs), 
      dK(n_regs, n_regs),
      d2K(n_regs, n_regs),
      opt_fix_effect(opt_fix_effect)
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
    sampleW_VY();
    }

    void sampleW_VY();
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
    }

// scale is too big, why
    double grad_theta_sigma_eps() const {
        double g = 0;
        if (family=="normal") {
            VectorXd tmp = Y - A * getW();
            g = log(sigma_eps) - pow(sigma_eps, -3) * tmp.dot(tmp) / n_obs;
        }
// std::cout << "n_obs=" << n_obs << std::endl;
        // return -g * sigma_eps / (n_obs*n_obs);
        return 0;
    }

    void set_theta_sgima_eps(double theta) {
        sigma_eps = exp(theta);
    }

// fixed effects
    VectorXd grad_beta() const {
        VectorXd grads = pow(sigma_eps, -2) * X.transpose() * (Y - X*beta - A*getW());
        return -grads / (n_obs * n_obs);
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
    VectorXd avg_gradient (n_paras);
        avg_gradient = VectorXd::Constant(n_paras, 0);
    
    for (int i=0; i < n_gibbs; i++) {
        
        // stack grad
        VectorXd gradient (n_paras); 
            gradient = VectorXd::Constant(n_paras, 0);
        
        int pos = 0;
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int theta_len = (*it)->getThetaSize();
auto timer_computeg = std::chrono::steady_clock::now();
            gradient.segment(pos, theta_len) = (*it)->getGrad();
std::cout << "timer_computeg (ms): " << since(timer_computeg).count() << std::endl;   
            pos += theta_len;
        }

        avg_gradient += gradient;
        sampleV_WY(); 
auto timer_sampleW = std::chrono::steady_clock::now();
        sampleW_VY();
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