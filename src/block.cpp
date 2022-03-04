#include "block.h"
#include <cmath>

using std::pow;

// ---- helper function for sampleW ----
Eigen::VectorXd rnorm_vec(int n, double mu, double sigma)
{
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; i++)
  {
    out[i] = R::rnorm(mu, sigma);
  }
  return (out);
}

// ---- other functions ------

void BlockModel::setW(const VectorXd& W) {
  int pos = 0;
  for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      int size = (*it)->getSize();
      (*it)->setW(W.segment(pos, size));
      pos += size;
  }
}

// sample W|VY 
void BlockModel::sampleW_VY()
{
  // for loop latent to build V
// mean <- 
// V <- sigma * V
  VectorXd SV = getSV();
  VectorXd inv_SV = VectorXd::Constant(SV.size(), 1).cwiseQuotient(SV);

  SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
  SparseMatrix<double> QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;

  chol_Q.compute(QQ);
  VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() + pow(sigma_eps,2) * A.transpose() * Y;

  VectorXd z (n_regs); 
  z = rnorm_vec(n_regs, 0, 1);
  
  // sample W ~ N(QQ^-1*M, QQ^-1)
  VectorXd W = chol_Q.rMVN(M, z);

// std::cout<< "V = " << V << std::endl;
// std::cout<< "W = " << W << std::endl;
  setW(W);
}

// std::cout << "W=" << W << std::endl;

// precondioner

  // hess_alpha <- function(alpha, V,  W){
  // n = n_obs
  // K_matrix <- diag(n)
  // K_matrix[seq(2, n*n, by=n+1)] <- -alpha
  // dK <-matrix(data=0, nrow=n, ncol=n)
  // dK[seq(2, n*n, by=n+1)] <- -1

  // M <- solve(as.matrix(K_matrix),dK)

  // Q <- t(K_matrix)%*%diag(1/V)%*%K_matrix
  // m <- solve(K_matrix, -1 + V)
  // Q_tilde <- Q + Asq
  // m_tilde <- solve(Q_tilde, Q%*%m + t(A_matrix)%*%Y)

  // dQ <- 1/sigma * (t(dK) %*%diag(1/V)%*%K_matrix +
  //                     t(K_matrix)%*%diag(1/V)%*%dK)

  // temp <- dQ%*%m - Q%*%forwardsolve(K_matrix,dK, upper.tri = F) %*% forwardsolve(K_matrix,-1 + V,upper.tri = F)

  // Q_t_inv <- solve(Q_tilde)

  // dm_tilde <- -Q_t_inv%*%dQ%*%
  //   Q_t_inv%*%(Q%*%m + sigma_eps^(-2)*t(A_matrix)%*%Y) +
  //   Q_t_inv%*%temp

  // term1 <- -sum(diag(M%*%M))
  // term2 <- sum(diag(t(dK)%*%diag(1/V)%*%dK%*%Q_t_inv)) -
  //   sum(diag(t(dK)%*%diag(1/V)%*%K_matrix%*%Q_t_inv%*%dQ%*%Q_t_inv))
  // term3 <- t(dm_tilde)%*%t(dK)%*%diag(1/V)%*%K_matrix%*%m_tilde +
  //   t(m_tilde)%*%t(dK)%*%diag(1/V)%*%dK%*%m_tilde +
  //   t(m_tilde)%*%t(dK)%*%diag(1/V)%*%K_matrix%*%dm_tilde
  // term4 <- t(dm_tilde)%*%t(dK)%*%diag(1/V)%*%(1-V)

  // return(term1-term2-term3-term4)


// }


