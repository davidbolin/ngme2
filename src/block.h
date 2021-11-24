#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <iostream>

#include "model.h"
#include "latent.h"

using Eigen::SparseMatrix;

class Block : public Model {
protected:
    unsigned n_obs, n_latent;

    std::vector<Latent> latents;

    VectorXd Y, W, V; 
    MatrixXd As, Ks, dKs;

    SparseMatrix<double> A, K, dK;
    // SparseMatrix<double> As, Ks;
public:
    Block() {}
    Block(VectorXd, Rcpp::List);
    ~Block() {}

    void sampleW(){};
    MatrixXd& precond();
    VectorXd& grad() {};
    VectorXd& get_parameter();
    void      set_parameter(const VectorXd&);

// //  m <- sparseSolve(K_matrix, -1 + V, upper.tri =F) chloskey
// // (sum(diag(solve(as.matrix(K_matrix),dK))) - 

    // horizontal
    void assembleA() {
        if (n_latent == 1) {
            As = latents[0].getA();
        }
    };  
    
    // block diagnol
    void assembleK() {
        // get K and assemble
        if (n_latent == 1) {
            Ks = latents[0].getK();
        }
    }  
    
    void assemble_dK() {
        if (n_latent == 1) {
            dKs = latents[0].get_dK();
        }
    }  

    // huge vector
    void assembleV() {
        if (n_latent == 1) {
            V = latents[0].getV();
        }
    }

    Rcpp::List testResult() {
        Rcpp::List res;

        res["A"] = As;
        res["K"] = Ks;
        res["dK"] = dKs;
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
Block::Block(VectorXd Y,
             Rcpp::List latents_in) 
    //init
    : n_obs(Y.size()), n_latent(latents_in.size()), Y(Y), 
      As(n_obs, n_latent * n_obs), Ks(n_latent * n_obs, n_latent * n_obs)
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
Block::_grad() {
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
Block::get_parameter() {
    for (std::vector<Latent>::iterator it = latents.begin(); it != latents.end(); it++) {
// concat theta_K, theta_m, theta_V
// stack together?
    }
}

inline void
Block::set_parameter(const VectorXd&) {
    for (std::vector<Latent>::iterator it = latents.begin(); it != latents.end(); it++) {
// set theta_K, theta_m, theta_V
    }
}

inline MatrixXd&
Block::precond() {
    // use previous W and V or new one

  MatrixXd a(1,1); return a;
}


#endif