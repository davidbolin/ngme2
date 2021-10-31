#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include <Eigen/Sparse>
#include <string>
#include <vector>
#include "latent.h"

using std::string;
using Eigen::SparseMatrix;

class Block : Model
{
friend class Latent;

private:
    VectorXd Y, w; 

    // string noise;
    Var noise;

    SparseMatrix<double> As, Ks;

    // std::vector<Model> models;
    std::vector<Latent> latents;
    
public:
    Block(){}
    Block(std::vector<Latent>) {} // -> construct As and Ks

    // to-do: generate new samples using current w
    void sample_w();


    // to-do: how to generate preconditioner
    MatrixXd& const precond();

    // to-do: generate gradient
    VectorXd& const grad();
    
    // to-do: build As into block.
    void assembleA(); 
    void assembleK();

    void set_theta(const VectorXd& theta); 
    // { dispatch theta according to latents }
};


inline VectorXd& const
Block::grad() {
    // call latent model gradient
    
    // vector<VectorXd> -> VectorXd

    VectorXd a (1);
    return a;
}



#endif