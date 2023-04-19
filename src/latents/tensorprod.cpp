// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x

// group="iid"
// group=ar1()

#include "../latent.h"
#include "MatrixAlgebra.h"

// tensor product for the C G class
Tensor_prod::Tensor_prod(const Rcpp::List& model_list, unsigned long seed)
  : Latent(model_list, seed) {
  // std::cout << "start init" << std::endl;

  Rcpp::List group = model_list["left"];
  left = LatentFactory::create(group, latent_rng());
  Rcpp::List model_right = model_list["right"];
  right = LatentFactory::create(model_right, latent_rng());
  n_theta_l = left->get_n_theta_K();
  n_theta_r = right->get_n_theta_K();
  V_size_l = left->get_V_size();
  V_size_r = right->get_V_size();
  W_size_l = left->get_W_size();
  W_size_r = right->get_W_size();

  K = getK(theta_K);
  for (int i=0; i < n_rep; i++)
    setSparseBlock(&K_rep, i*V_size, i*W_size, K);

  SparseMatrix<double> Q = K.transpose() * K;
  if (W_size == V_size && left->symmetricK && right->symmetricK) {
    symmetricK = true;
    chol_solver_K.init(W_size, 0,0,0);
    chol_solver_K.analyze(K);
  } else {
    lu_solver_K.init(W_size, 0,0,0);
    lu_solver_K.analyze(K);
  }

  // Init QV
  solver_Q.init(W_size, 0,0,0);
  solver_Q.analyze(Q);
  update_each_iter();

  // std::cout << "Finish init" << std::endl;
}

SparseMatrix<double> Tensor_prod::getK(const VectorXd& theta_K) const {
  SparseMatrix<double> Kl = left->getK(theta_K.segment(0, n_theta_l));
  SparseMatrix<double> Kr = right->getK(theta_K.segment(n_theta_l, n_theta_r));

  return kronecker(Kl, Kr);
}

SparseMatrix<double> Tensor_prod::get_dK(int index, const VectorXd& alpha) const {
  VectorXd theta_K_l = theta_K.segment(0, n_theta_l);
  VectorXd theta_K_r = theta_K.segment(n_theta_l, n_theta_r);

  if (index < n_theta_l) {
    SparseMatrix<double> dk_l = left->get_dK(index, theta_K_l);
    SparseMatrix<double> K_r = right->getK(theta_K_r);
    return kronecker(dk_l, K_r);
  } else {
    SparseMatrix<double> dk_r = right->get_dK(index - n_theta_r, theta_K_r);
    SparseMatrix<double> K_l = left->getK(theta_K_l);
    return kronecker(K_l, dk_r);
  }
}

// ------ iid model -----

Iid::Iid(const Rcpp::List& model_list, unsigned long seed)
  : Latent(model_list, seed),
    I(Rcpp::as< SparseMatrix<double,0,int> > (model_list["K"])) {
    K = I;
    fix_flag[latent_fix_theta_K] = true;
}

SparseMatrix<double> Iid::getK(const VectorXd& alpha) const {return I;};
SparseMatrix<double> Iid::get_dK(int index, const VectorXd& alpha) const {return 0*I;};

