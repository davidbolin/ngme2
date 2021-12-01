#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <iostream>

#include "include/MatrixAlgebra.h"
#include "model.h"
#include "latent.h"

using Eigen::SparseMatrix;

class BlockModel : public Model {
protected:
    unsigned n_obs, n_latent;
    std::vector<Latent> latents;

    VectorXd Y, W, V, Mu;
    SparseMatrix<double> A, K, dK, d2K;
public:
    BlockModel() {}
    BlockModel(VectorXd, Rcpp::List);
    ~BlockModel() {}

    void sampleW(){};
    MatrixXd& precond();
    VectorXd& grad() {};
    VectorXd& get_parameter();
    void      set_parameter(const VectorXd&);

// //  m <- sparseSolve(K_matrix, -1 + V, upper.tri =F) chloskey
// // (sum(diag(solve(as.matrix(K_matrix),dK))) - 

    // horizontal
    void assembleA() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&A, 0, i*n_obs, latents[i].getA());
        }
    };  
    
    // block diagnol
    void assembleK() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&K, i*n_obs, i*n_obs, latents[i].getK());
        }
    }  
    
    void assemble_dK() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&dK, i*n_obs, i*n_obs, latents[i].get_dK());
        }
    }  

    // huge vector
    void assembleV() {
        // if (n_latent == 1) {
        //     V = latents[0].getV();
        // }
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&K, i*n_obs, i*n_obs, latents[i].get_dK());
        }
    }
    
    // huge vector
    void assembleW() {
        for (unsigned i=0; i < n_latent; i++) {
            W.segment(i*n_obs, i*(n_obs+1)) = latents[i].getW();
        }
    }

    // huge vector
    void assembleMu() {
        for (unsigned i=0; i < n_latent; i++) {
            Mu.segment(i*n_obs, i*(n_obs+1)) = latents[i].getMu();
        }
    }

    void assembleGrad() {

    }    

    Rcpp::List testResult() {
        Rcpp::List res;

        res["A"] = A;
        res["K"] = K;
        res["dK"] = dK;
        res["V"] = V;
        return res;
    }

//     // void getA(Latent l) {};
//     // void set_theta(const VectorXd& theta); 
//     // { dispatch theta according to latents }


    VectorXd& _grad();
//     // VectorXd& _grad_rb();
};


inline 
BlockModel::BlockModel(VectorXd Y,
             Rcpp::List latents_in) 
    //init
    : n_obs(Y.size()), n_latent(latents_in.size()), Y(Y), 
      A(n_obs, n_latent * n_obs), K(n_latent * n_obs, n_latent * n_obs), V(n_latent * n_obs)
    {
// std::cout << latents_in.size() << std::endl;

    for (int i=0; i < n_latent; ++i) {
        Rcpp::List latent_in = latents_in[i];

        // construct acoording to different models
        string type = latent_in["type"];
        if (type == "ar1") {
            latents.push_back( AR(latent_in) );
        }
    }
    
    assembleA();
    assembleK();
    assembleV();
    assemble_dK();

    // assenble_mu();
}

inline VectorXd&
BlockModel::_grad() {
    // solve(K) sparse
    // chloskey solver for sparse matrix

    Eigen::SparseLU<SparseMatrix<double> > solver; 
    solver.analyzePattern(K);   // for this step the numerical values of A are not used
    solver.factorize(K);
    
    // solver.solve(dK).trace(); // solve a matrix
    // -
    // W.transpose() * (dK) * (VectorXd::Constant(V.size(), 1).cwiseQuotient(V).asDiagonal() * (K * W + (VectorXd::Constant(V.size(), 1) - V) * mu);
    
    //// (sum(diag(solve(K, dK))) - 
    ////    t(W) %*% dK %*% diag(1/V) %*%(K_matrix%*%W+(1-V) * mu))

    // VectorXd a (1); return a;
}

inline VectorXd& 
BlockModel::get_parameter() {
    for (std::vector<Latent>::iterator it = latents.begin(); it != latents.end(); it++) {
// concat theta_K, theta_m, theta_V
// stack together?
    }
}

inline void
BlockModel::set_parameter(const VectorXd&) {
    for (std::vector<Latent>::iterator it = latents.begin(); it != latents.end(); it++) {
// set theta_K, theta_m, theta_V
    }
}

inline MatrixXd&
BlockModel::precond() {
    // use previous W and V or new one

  MatrixXd a(1,1); return a;
}


#endif