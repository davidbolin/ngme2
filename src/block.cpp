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
// if (debug) std::cout << "starting sampling W." << std::endl;        

  VectorXd SV = getSV();
  VectorXd inv_SV = VectorXd::Constant(SV.size(), 1).cwiseQuotient(SV);

  SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
  SparseMatrix<double> QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;

  chol_QQ.compute(QQ);
  VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() + 
      pow(sigma_eps, -2) * A.transpose() * (Y - X * beta);

  VectorXd z (n_meshs); 
  z = rnorm_vec(n_meshs, 0, 1);
  
  // sample W ~ N(QQ^-1*M, QQ^-1)
  VectorXd W = chol_QQ.rMVN(M, z);
// if (debug) std::cout << "after rMVN" << std::endl;        

  setW(W);

// fixW
if (fixW) setW(trueW);
  
// if (debug) std::cout << "finish sampling W." << std::endl;        
}

// sample W|V
void BlockModel::sampleW_V()
{
  // sample KW ~ N(mu*(V-h), diag(V))
  VectorXd SV = getSV();
  Eigen::VectorXd KW (n_meshs);
  for (int i=0; i < n_meshs; i++) {
    KW[i] = R::rnorm(0, sqrt(SV[i]));
  }
  KW = getMean() + KW;

  LU_K.factorize(K);
  VectorXd W = LU_K.solve(KW);

  setW(W);
}


