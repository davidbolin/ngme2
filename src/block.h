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

    VectorXd Y, w, V; 
    MatrixXd As, Ks;

    // SparseMatrix<double> As, Ks;
    
public:
    Block() {}
    Block(VectorXd, Rcpp::List);
    ~Block() {}

//     // void sample_w();
    MatrixXd& precond();
    VectorXd& grad();
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
    
    // huge vector
    void assembleV() {
        if (n_latent == 1) {
            V = latents[0].getV();
        }
    }

    Rcpp::List result() {
        Rcpp::List res;

        res["A"] = As;
        res["K"] = Ks;
        res["V"] = V;
        return res;
    }


//     // void getA(Latent l) {};

//     // void set_theta(const VectorXd& theta); 
//     // { dispatch theta according to latents }


//     // VectorXd& _grad();
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
// std::cout << "yes" << std::endl;
            latents.push_back( AR(latent_in) );
        } else {
            
        }
    }
    
    assembleA();
    assembleK();
    assembleV();
}

inline VectorXd&
Block::grad() {
    // use current w and V
    

    VectorXd a (1); return a;
}

inline MatrixXd&
Block::precond() {
    // use previous w and V or new one

  MatrixXd a(1,1); return a;
}

inline VectorXd& 
Block::get_parameter() {
    VectorXd a (1); return a;
}

inline void
Block::set_parameter(const VectorXd&) {

}

#endif