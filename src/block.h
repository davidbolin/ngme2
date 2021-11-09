#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include <Eigen/Sparse>
#include <string>
#include <vector>
#include "model.h"
#include "latent.h"

using std::string;
using Eigen::SparseMatrix;

class Block : public Model
{
friend class Latent;
protected:
    VectorXd Y, w; 
    SparseMatrix<double> As, Ks;

    string noise; // Var noise;

    std::vector<Latent> latents;
    
public:
    Block(){}
    Block(std::vector<Latent>) {} // -> construct As and Ks

    void sample_w();

    MatrixXd& const precond();
    VectorXd& const grad(); // call _grad() or _grad_rb()

//  m <- sparseSolve(K_matrix, -1 + V, upper.tri =F) chloskey
// (sum(diag(solve(as.matrix(K_matrix),dK))) - 

    // to-do: build block
    void assembleA();  // horizontal
    void assembleK();  // block diagnol
    void aseembleV();  // huge vector

    void getA(Latent l) {};

    void set_theta(const VectorXd& theta); 
    // { dispatch theta according to latents }

    VectorXd& const _grad();
    VectorXd& const _grad_rb();
};


inline VectorXd& const
Block::_grad() {
    // (sum(diag(solve(as.matrix(K_matrix), dK))) - 
    //    t(W) %*% dK %*% diag(1/V) %*%(K_matrix%*%W+(1-V) * mu))

    VectorXd a (1); return a;
}

inline VectorXd& const
Block::_grad_rb() {
    // sum(diag(solve(Q_tilde,t(dK)%*%diag(1/V)%*%K_matrix)))-
    // t(m_tilde)%*%t(dK) %*%diag(1/V)%*%K_matrix%*%m_tilde - 
    // t(m_tilde)%*%t(dK)%*%diag(1/V)%*%(1-V)*mu)

    VectorXd a (1); return a;
}

void sampleW(){ 
    // A_obs, huge A, K <- huge K , V <- huge V
    Eigen::SparseMatrix<double> Asq = A_obs.transpose() * A_obs / (sigma_eps * sigma_eps);


    Eigen::SparseMatrix<double> Q = K * inv_V.asDiagonal() * K; 

    Q = Q/(sigma_eps*sigma_eps);

    Q = Q + Asq; 

    Eigen::VectorXd resp = K * inv_V.asDiagonal() * (-h + V_temp);

    resp = resp.cwiseProduct(mu_gibbs)/(sigma_eps*sigma_eps);

    resp = resp + A_obs.transpose()*Y/(sigma_eps * sigma_eps);

// chloskey solver for sparse matrix
    Eigen::SimplicialLLT <Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> chol_Q;
    chol_Q.analyzePattern(Q);
    chol_Q.factorize(Q);

    Eigen::VectorXd m_W_new = chol_Q.solve(resp);
//

    Eigen::MatrixXd chol_W;

    chol_W = chol_Q.matrixU().eval();

    m_W.col(i+1) = m_W_new;

    Eigen::VectorXd norm_temp;

    norm_temp = rnorm_vec(m,0,1);

    Eigen::VectorXd new_W;

    new_W = chol_W.triangularView<Eigen::Upper>().solve(norm_temp);

    new_W = m_W_new + new_W;
}

MatrixXd& const
precond() {
  hess_alpha <- function(alpha, V,  W){
  n = n_obs
  K_matrix <- diag(n)
  K_matrix[seq(2, n*n, by=n+1)] <- -alpha
  dK <-matrix(data=0, nrow=n, ncol=n)
  dK[seq(2, n*n, by=n+1)] <- -1
  
  M <- solve(as.matrix(K_matrix),dK)
  
  Q <- t(K_matrix)%*%diag(1/V)%*%K_matrix
  m <- solve(K_matrix, -1 + V)
  Q_tilde <- Q + Asq
  m_tilde <- solve(Q_tilde, Q%*%m + t(A_matrix)%*%Y)
  
  dQ <- 1/sigma * (t(dK) %*%diag(1/V)%*%K_matrix + 
                      t(K_matrix)%*%diag(1/V)%*%dK)
  
  temp <- dQ%*%m - Q%*%forwardsolve(K_matrix,dK, upper.tri = F) %*% forwardsolve(K_matrix,-1 + V,upper.tri = F)
  
  Q_t_inv <- solve(Q_tilde)
  
  dm_tilde <- -Q_t_inv%*%dQ%*%
    Q_t_inv%*%(Q%*%m + sigma_eps^(-2)*t(A_matrix)%*%Y) +
    Q_t_inv%*%temp
  
  term1 <- -sum(diag(M%*%M))
  term2 <- sum(diag(t(dK)%*%diag(1/V)%*%dK%*%Q_t_inv)) -
    sum(diag(t(dK)%*%diag(1/V)%*%K_matrix%*%Q_t_inv%*%dQ%*%Q_t_inv))
  term3 <- t(dm_tilde)%*%t(dK)%*%diag(1/V)%*%K_matrix%*%m_tilde +
    t(m_tilde)%*%t(dK)%*%diag(1/V)%*%dK%*%m_tilde +
    t(m_tilde)%*%t(dK)%*%diag(1/V)%*%K_matrix%*%dm_tilde
  term4 <- t(dm_tilde)%*%t(dK)%*%diag(1/V)%*%(1-V)
  
  return(term1-term2-term3-term4)
}
}

#endif