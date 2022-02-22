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

#include "include/solver.h"
#include "include/MatrixAlgebra.h"
#include "model.h"
#include "latent.h"

#include "latents/ar1.h"

using Eigen::SparseMatrix;

class BlockModel : public Model {
protected:
// n_regs   = row(A1) + ... + row(An)
// n_paras = total paras need to optimize
    unsigned n_obs, n_latent, n_regs, n_paras;
    std::vector<Latent*> latents;

// Q: what is mu, sigma here
    double mu, sigma, sigma_eps;
    
    VectorXd Y, Mu, h;

    SparseMatrix<double> A, K, dK, d2K;

    cholesky_solver chol_Q;
public:
    BlockModel() {}
    // main constructor
    BlockModel(const VectorXd Y, 
               const int n_paras,
               const int n_regs,
               Rcpp::List latents_in) 
    : n_obs(Y.size()),
      n_latent(latents_in.size()), 
      n_regs(n_regs),
      n_paras(n_paras),
      mu(1),
      sigma(1),
      sigma_eps(1),
      h(n_regs),
      Y(Y), 
      A(n_obs, n_regs), 
      K(n_regs, n_regs), 
      dK(n_regs, n_regs),
      d2K(n_regs, n_regs),
      Mu(n_regs)
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
    
    /* Init variables */
    h = VectorXd::Constant(n_regs, 1);

    assemble();
        VectorXd inv_V = VectorXd::Constant(n_regs, 1).cwiseQuotient(getV());
        SparseMatrix<double> QQ = pow(sigma, -2) * K.transpose() * inv_V.asDiagonal() * K + pow(sigma_eps, 2) * A.transpose() * A;
    chol_Q.analyze(QQ);
}

    ~BlockModel() {}

    void sampleW_VY();
    void sampleV();
    void setW(const VectorXd&);

    // void gibbsSample();

    VectorXd             get_parameter() const;
    void                 set_parameter(const VectorXd&);
    VectorXd             grad() const;
    SparseMatrix<double> precond() const;

    void assemble() {
        int n = 0;
        for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
            setSparseBlock(&A,   0, n, (*it)->getA());         
            setSparseBlock(&K,   n, n, (*it)->getK());      
            setSparseBlock(&dK,  n, n, (*it)->get_dK());   
            setSparseBlock(&d2K, n, n, (*it)->get_d2K()); 
            
            n += (*it)->getSize();
        }
    }

    VectorXd getV() const {
        VectorXd V (n_regs);
        int pos = 0;
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int size = (*it)->getSize();
            V.segment(pos, size) = (*it)->getV();
            pos += size;
        }
        return V;
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

    // void assembleTheta() {
    //     int pos = 0;
    //     for (unsigned i=0; i < n_latent; i++) {
    //         int theta_len = (*latents[i]).getTheta().size();
    //         Theta.segment(pos, theta_len) = (*latents[i]).getTheta();
    //         pos += theta_len;
    //     }
    // }

    // void assembleGrad() {
    //     int pos = 0;
    //     for (unsigned i=0; i < n_latent; i++) {
    //         int theta_len = (*latents[i]).getTheta().size();
    //         Grad.segment(pos, theta_len) = (*latents[i]).getGrad();
    //         pos += theta_len;
    //     }
    // }    

//     // void getA(Latent l) {};
//     // void set_theta(const VectorXd& theta); 
//     // { dispatch theta according to latents }

    // VectorXd& _grad_rb();


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
    res["V"] = getV();
    res["W"] = getW();
    return res;
}

void setVW() {
    // if (n_regs==10) {
    //     V << 0.9375010, 0.2752556, 0.3895697, 1.7309434, 0.5218133, 0.3935415, 1.7542207, 0.5890800, 0.6039723, 2.6892251;
    //     W << -0.63479296, -1.27376433, -1.35513971, -1.55145669, -0.87165923, -0.31881076, 1.04077067, 1.72322077,  0.01019182, 4.14311608;
    // }
};

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
    return thetas;
}

inline void BlockModel::set_parameter(const VectorXd& Theta) {
std::cout << "Theta=" << Theta <<std::endl;
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        int theta_len = (*it)->getThetaSize();
        VectorXd theta = Theta.segment(pos, theta_len);
        (*it)->setTheta(theta);
        pos += theta_len;
    }
}

inline VectorXd BlockModel::grad() const {
    VectorXd gradient (n_paras);

    return gradient;
}

inline SparseMatrix<double> BlockModel::precond() const {
    SparseMatrix<double> precond;

    return precond;
}



#endif