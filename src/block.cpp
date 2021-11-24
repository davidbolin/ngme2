#include "block.h"

inline void 
Block::sampleW(){ 
    // A, huge A, K <- huge K , V <- huge V
    VectorXd inv_V = VectorXd::Constant(V.size(), 1).cwiseQuotient(V);
    
    Eigen::SparseMatrix<double> Asq = A.transpose() * A / (sigma_eps * sigma_eps);
    Eigen::SparseMatrix<double> Q = K * inv_V.asDiagonal() * K; 

    Q = Q/(sigma_eps*sigma_eps);

    Q = Q + Asq; 

    Eigen::VectorXd resp = K * inv_V.asDiagonal() * (-h + V_temp);

    resp = resp.cwiseProduct(mu_gibbs)/(sigma_eps*sigma_eps);

    resp = resp + A.transpose()*Y/(sigma_eps * sigma_eps);

// chloskey solver for sparse matrix
    Eigen::SimplicialLLT <Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > chol_Q;
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

MatrixXd&
precond() {
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

  MatrixXd a(1,1); return a;
}
