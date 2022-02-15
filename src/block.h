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
    unsigned n_obs, n_latent;
    std::vector<Latent*> latents;
    double mu=1;

    SparseLU<SparseMatrix<double> > solver;
    VectorXd Y, W, V, Mu, Grad, Theta;
    SparseMatrix<double> A, K, dK, d2K;
    
    // lu_sparse_solver K_solver;
public:
    BlockModel() {}
    BlockModel(VectorXd, Rcpp::List);
    ~BlockModel() {}

    void sampleW();
    void setW(const VectorXd&);

    void gibbsSample();

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


    void _grad();
    // VectorXd& _grad_rb();

    Rcpp::List testResult();
    Rcpp::List testGrad();
};

// ------------------ Constructor ----------------------

inline 
BlockModel::BlockModel(VectorXd Y, Rcpp::List latents_in) 
    : n_obs(Y.size()), 
      n_latent(latents_in.size()), Y(Y), 
      A(n_obs, n_latent * n_obs), 
      K(n_latent * n_obs, n_latent * n_obs), 
      dK(n_latent * n_obs, n_latent * n_obs),
      d2K(n_latent * n_obs, n_latent * n_obs),
      V(n_latent * n_obs), 
      W(n_latent * n_obs),
      Mu(n_latent * n_obs),
      Theta(n_latent), Grad(n_latent)
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

    solver.analyzePattern(K);
    solver.factorize(K);
}




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