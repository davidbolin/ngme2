#include "block.h"

// helper function for sampleW
Eigen::VectorXd rnorm_vec(int n, double mu, double sigma)
{
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; i++)
  {
    out[i] = R::rnorm(mu, sigma);
  }
  return (out);
}


inline void
BlockModel::set_parameter(const VectorXd& Theta) {
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        int theta_len = (*it)->getTheta().size();
        VectorXd theta = Theta.segment(pos, pos + theta_len);
// std::cout << theta << std::endl;
        (*it)->setTheta(theta);
        pos += theta_len;
    }
}


inline VectorXd &
BlockModel::grad() {
  assembleK();
  _grad();
  return Grad;
}

// compute gradient 
inline void
BlockModel::_grad()
{
  VectorXd h = VectorXd::Constant(V.size(), 1);
  
  // 1. To compute tr(dK * K^-1)
  SparseMatrix<double> I (n_obs, n_obs); 
  I.setIdentity();
  double g = (dK * solver.solve(I)).eval().diagonal().sum();

  // 2. Compute the rest
  double rhs = W.transpose() * dK.transpose() * 
               (VectorXd::Constant(n_obs, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * mu);
  
// std::cout << g-rhs << std::endl;
  Grad << g - rhs;
}


// sample W|V
inline void
BlockModel::sampleW()
{
  std::cout << "V =" << V << std::endl;
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

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> chol_Q;
  chol_Q.analyzePattern(Q);
  chol_Q.factorize(Q);

  Eigen::VectorXd m_W_new = chol_Q.solve(resp);

  Eigen::MatrixXd chol_W = chol_Q.matrixU().eval();

  Eigen::VectorXd norm_temp = rnorm_vec(V_len, 0, 1);

  Eigen::VectorXd new_W = chol_W.triangularView<Eigen::Upper>().solve(norm_temp);

  new_W = m_W_new + new_W;

  // distribute W
  setW(new_W);

  // std::cout << "W =" << new_W << std::endl;
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

  sampleW();

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
  res["W"] = W;

  
  // res["Mu"] = Mu;
  // ind_IG ig(n_obs, 0.5, 0.5);
  // ig.sample_cond_V((*latents[0]).getW(), Mu, K);
  
  // grad();
// std::cout << Grad << std::endl;

  return res;
}


Rcpp::List
BlockModel::testGrad()
{
  // set TrueV
  // std::cout << "V=" << V.size() << std::endl;
  // std::cout << "W=" << W.size() << std::endl;
  V << 0.9375010, 0.2752556, 0.3895697, 1.7309434, 0.5218133, 0.3935415, 1.7542207, 0.5890800, 0.6039723, 2.6892251;
  W << -0.63479296, -1.27376433, -1.35513971, -1.55145669, -0.87165923, -0.31881076, 1.04077067, 1.72322077,  0.01019182, 4.14311608;

  Rcpp::List res;
  return res;
}

// std::cout << "Q =" << Q << std::endl;
// std::cout << "K =" << K << std::endl;
// std::cout << "resp1 =" << resp << std::endl;
// std::cout << "resp2 =" << resp << std::endl;