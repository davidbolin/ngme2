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

// 

class BlockModel : public Model {
protected:
    unsigned n_obs, n_latent;
    std::vector<Latent*> latents;

    VectorXd Y, W, V, Mu, Grad, Theta;
    SparseMatrix<double> A, K, dK, d2K;
    
    lu_sparse_solver K_solver;
public:
    BlockModel() {}
    BlockModel(VectorXd, Rcpp::List);
    ~BlockModel() {}

    // W|V
    void sampleW();
    void gibbsSample();

    VectorXd& get_parameter();
    void set_parameter(const VectorXd&);
    MatrixXd& precond();
    VectorXd& grad();

// //  m <- sparseSolve(K_matrix, -1 + V, upper.tri =F) chloskey
// // (sum(diag(solve(as.matrix(K_matrix),dK))) - 

    // horizontal
    void assembleA() {
// setSparseBlock(&A, 0, 0*n_obs, latents[0].getA());
// std::cout << A << std::endl;
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&A, 0, i*n_obs, (*latents[i]).getA());
        }
    };  
    
    // block diagonal
    void assembleK() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&K, i*n_obs, i*n_obs, (*latents[i]).getK());
        }
    }  
    
    void assemble_dK() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&dK, i*n_obs, i*n_obs, (*latents[i]).get_dK());
        }
    }

    void assemble_d2K() {
        for (unsigned i=0; i < n_latent; i++) {
            setSparseBlock(&d2K, i*n_obs, i*n_obs, (*latents[i]).get_d2K());
        }
    }  

    // vector
    void assembleV() {
        for (unsigned i=0; i < n_latent; i++) {
            V.segment(i*n_obs, n_obs) = (*latents[i]).getV();
        }
    }
    
    void assembleW() {
        for (unsigned i=0; i < n_latent; i++) {
            W.segment(i*n_obs, n_obs) = (*latents[i]).getW();
        }
    }

    void assembleMu() {
        for (unsigned i=0; i < n_latent; i++) {
            Mu.segment(i*n_obs, n_obs) = (*latents[i]).getMu();
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


    VectorXd& _grad();
    // VectorXd& _grad_rb();

    Rcpp::List testResult();
};

// ---- Constructor ----

inline 
BlockModel::BlockModel(VectorXd Y, Rcpp::List latents_in) 
    : n_obs(Y.size()), 
      n_latent(latents_in.size()), Y(Y), 
      A(n_obs, n_latent * n_obs), 
      K(n_latent * n_obs, n_latent * n_obs), 
      dK(n_latent * n_obs, n_latent * n_obs),
      d2K(n_latent * n_obs, n_latent * n_obs),
      V(n_latent * n_obs), 
      Mu(n_latent * n_obs),
      Theta(30), Grad(30)
    {

    for (int i=0; i < n_latent; ++i) {
        Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);

        // construct acoording to different models
        string type = latent_in["type"];
        if (type == "ar1") {
            latents.push_back(new AR(latent_in) );
        }
    }
    
    Mu = VectorXd::Zero(n_latent * n_obs);
    assembleA();
    assembleK();
    assemble_dK();
    assemble_d2K();
    assembleV();
    assembleMu();
}

// ---- inherited functions ------

inline MatrixXd&
BlockModel::precond() {
    MatrixXd a(1,1); return a;
}

inline VectorXd& 
BlockModel::get_parameter() {
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        VectorXd theta = (*it)->getTheta();
        Theta.segment(pos, pos + theta.size()) = theta;
        pos += theta.size();
    }
    return Theta;
}

inline VectorXd& 
BlockModel::grad() {
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        VectorXd grad = (*it)->getGrad();
        Grad.segment(pos, pos + grad.size()) = grad;
        pos += grad.size();
    }
    return Grad;
}

inline void
BlockModel::set_parameter(const VectorXd&) {
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        int theta_len = (*it)->getTheta().size();
        VectorXd theta = Theta.segment(pos, pos + theta_len);
        (*it)->setTheta(theta);
        pos += theta_len;
    }
}


#endif