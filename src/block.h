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
// n_reg   = row(A1) + ... + row(An)
// n_paras = total paras need to optimize
    unsigned n_obs, n_latent, n_reg, n_paras;
    std::vector<Latent*> latents;

    double mu;
    VectorXd Y, W, V, Mu, Grad, Theta;
    SparseMatrix<double> A, K, dK, d2K;

    SparseLU<SparseMatrix<double> > solver;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> chol_Q;
    
    // lu_sparse_solver K_solver;
public:
    BlockModel() {}
    // main constructor
    BlockModel(const VectorXd Y, 
               const int n_paras,
               const int n_reg,
               Rcpp::List latents_in) 
    : n_obs(Y.size()),
      n_latent(latents_in.size()), 
      n_reg(n_reg),
      n_paras(n_paras),
      Y(Y), 
      A(n_obs, n_reg), 
      K(n_reg, n_reg), 
      dK(n_reg, n_reg),
      d2K(n_reg, n_reg),
      V(n_reg), 
      W(n_reg),
      Mu(n_reg),
      Theta(n_latent), 
      Grad(n_latent)
    {
    // Init each latent model
    for (int i=0; i < n_latent; ++i) {
        Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);

        // construct acoording to different models
        string type = latent_in["type"];
        if (type == "ar1") {
            latents.push_back(new AR(latent_in) );
        }
    }
    
    assemble();
    
    // VectorXd inv_V = VectorXd::Constant(V.size(), 1).cwiseQuotient(V);
    // Eigen::SparseMatrix<double> Q = K.transpose() * inv_V.asDiagonal() * K;
    // chol_Q.analyzePattern(Q);
    
    // solver.analyzePattern(K);
    // solver.factorize(K);
}

    ~BlockModel() {}

    void sampleW();
    void sampleV();
    void setW(const VectorXd&);

    // void gibbsSample();

    VectorXd& get_parameter();
    void set_parameter(const VectorXd&);
    MatrixXd& precond();
    VectorXd& grad();

    void assemble() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&A, 0, i*n_obs, (*latents[i]).getA());            // A
            setSparseBlock(&K, i*n_obs, i*n_obs, (*latents[i]).getK());      // K
            setSparseBlock(&dK, i*n_obs, i*n_obs, (*latents[i]).get_dK());   // dK
            setSparseBlock(&d2K, i*n_obs, i*n_obs, (*latents[i]).get_d2K()); // d2K
            // V.segment(i*n_obs, n_obs) = (*latents[i]).getV();
            // W.segment(i*n_obs, n_obs) = (*latents[i]).getW();
            // Mu.segment(i*n_obs, n_obs) = (*latents[i]).getMu();
        }
    }

    void assembleV() {
        for (unsigned i=0; i < n_latent; i++) {
            V.segment(i*n_obs, n_obs) = (*latents[i]).getV();
            // Mu.segment(i*n_obs, n_obs) = (*latents[i]).getMu();
        }
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


    double _grad();
    // VectorXd& _grad_rb();

    Rcpp::List testResult();
    Rcpp::List testGrad();
};

// inline VectorXd& 
// BlockModel::grad() {
//     int pos = 0;
//     for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
//         VectorXd grad = (*it)->getGrad();
//         Grad.segment(pos, pos + grad.size()) = grad;
//         pos += grad.size();
//     }
//     return Grad;
// }



#endif