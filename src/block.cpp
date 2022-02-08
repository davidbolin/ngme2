#include "block.h"
#include "include/solver.h"

Eigen::VectorXd rnorm_vec(int n, double mu, double sigma)
{
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; i++)
  {
    out[i] = R::rnorm(mu, sigma);
  }
  return (out);
}

inline VectorXd&
BlockModel::_grad() {
    // solve(K) sparse
    // chloskey solver for sparse matrix

    lu_sparse_solver K_solver;
    K_solver.analyze(K);
    K_solver.compute(K);

    VectorXd h = VectorXd::Constant(V.size(), 1);
    K_solver.trace(dK) - W.transpose() * dK.transpose() * (VectorXd::Constant(V.size(), 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * Mu);
    

    // Eigen::SparseLU<SparseMatrix<double> > solver; 
    // solver.analyzePattern(K);   // for this step the numerical values of A are not used
    // solver.factorize(K);
    
//Q: How to solve for a matrix? 
    // K_solver.solve(dK).trace(); // solve a matrix

    
    //// (sum(diag(solve(K, dK))) - 
    ////    t(W) %*% dK %*% diag(1/V) %*%(K_matrix%*%W+(1-V) * mu))

    // VectorXd a (1); return a;
}

// inline void
// Rcpp::List gibbs_sample(Eigen::VectorXd mu_coef, double kap,
//                         int N_sim, 
//                         Eigen::VectorXd V_init,
//                         Eigen::VectorXd W_init, Eigen::SparseMatrix<double> C,
//                         Eigen::SparseMatrix<double> G,
//                         Eigen::VectorXd h,
//                         Eigen::SparseMatrix<double> A_obs,
//                         Eigen::VectorXd Y,
//                         double eta,
//                         double sigma_eps, double sigma, Eigen::VectorXd loc
// ){
  
// }


// sample W|V
inline void
BlockModel::sampleW()
{
  VectorXd h = VectorXd::Constant(V.rows(), 1);

  double sigma_eps = 1;
  int V_len = V.rows();

  // A, K, V
  VectorXd inv_V = VectorXd::Constant(V.size(), 1).cwiseQuotient(V); 
  Eigen::SparseMatrix<double> Q = K.transpose() * inv_V.asDiagonal() * K;
  Eigen::SparseMatrix<double> Asq = A.transpose() * A / (sigma_eps * sigma_eps);

  Q = Q / (sigma_eps * sigma_eps);
  Q = Q + Asq;

  Eigen::VectorXd resp = K * inv_V.asDiagonal() * (-h + V);
  resp = resp.cwiseProduct(Mu) / (sigma_eps * sigma_eps);
  resp = resp + A.transpose() * Y / (sigma_eps * sigma_eps);


  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > chol_Q;
  chol_Q.analyzePattern(Q);
  chol_Q.factorize(Q);

  Eigen::VectorXd m_W_new = chol_Q.solve(resp);

  Eigen::MatrixXd chol_W;

  chol_W = chol_Q.matrixU().eval();

  Eigen::VectorXd norm_temp;
  norm_temp = rnorm_vec(V_len, 0, 1);

  Eigen::VectorXd new_W;
  new_W = chol_W.triangularView<Eigen::Upper>().solve(norm_temp);

  new_W = m_W_new + new_W;

std::cout << "W =" << new_W;
}

MatrixXd &
precond()
{
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

  MatrixXd a(1, 1);
  return a;
}

// for testing

Rcpp::List
BlockModel::testResult()
{

  Rcpp::List res;
A.makeCompressed();
K.makeCompressed();
dK.makeCompressed();
d2K.makeCompressed();
  res["A"] = A;
  res["K"] = K;
  res["dK"] = dK;
  res["d2K"] = d2K;
  res["V"] = V;
  res["Mu"] = Mu;

  sampleW();

  return res;
}