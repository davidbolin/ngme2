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
  std::cout << "start init" << std::endl;

  unsigned long latent_seed = latent_rng();
  Rcpp::List group = model_list["group"];
  left = LatentFactory::create(group, latent_seed);
  std::cout << "create left" << std::endl;
  Rcpp::List model_right = model_list["model_right"];
  right = LatentFactory::create(model_right, latent_seed);
  std::cout << "create right" << std::endl;

  K = getK(theta_K);
  for (int i=0; i < n_rep; i++)
    setSparseBlock(&K_rep, i*V_size, i*W_size, K);

  SparseMatrix<double> Q = K.transpose() * K;
  if (W_size == V_size) {
    lu_solver_K.init(W_size, 0,0,0);
    lu_solver_K.analyze(K);
  }

  // Init QV
  solver_Q.init(W_size, 0,0,0);
  solver_Q.analyze(Q);

  std::cout << "Finish init" << std::endl;
}

VectorXd Tensor_prod::grad_theta_K() {
  if (numer_grad)
    return numerical_grad();

  VectorXd grad = VectorXd::Zero(n_theta_K);
  int n1 = left->get_n_theta_K();
  int n2 = right->get_n_theta_K();

  SparseMatrix<double> K1 = left->getK();
  SparseMatrix<double> K2 = right->getK();

  // to-do later for non-tationary case
  SparseMatrix<double> dKl = left->get_dK_by_index(0);
  SparseMatrix<double> dKr = right->get_dK_by_index(0);
  SparseMatrix<double> dK  = kronecker(dKl, K2);
  SparseMatrix<double> dK2 = kronecker(K1, dKr);

  // compute trace manually
    double trace, trace2;
    if (!symmetricK) {
        lu_solver_K.computeKTK(K);
        trace = lu_solver_K.trace(dK);
        trace2 = lu_solver_K.trace(dK2);
    } else {
        chol_solver_K.compute(K);
        trace = chol_solver_K.trace(dK);
        trace2 = chol_solver_K.trace(dK2);
    }

  grad.segment(0, n1) = left->grad_theta_K(
    K, dK, Ws, prevWs, vars, mu, sigma, h, trace, W_size
  );

  grad.segment(n1, n2) = right->grad_theta_K(
    K, dK2, Ws, prevWs, vars, mu, sigma, h, trace2, W_size
  );

  return grad;
}

SparseMatrix<double> Tensor_prod::getK(const VectorXd& theta_K) const {
  int n1 = left->get_n_theta_K();
  int n2 = right->get_n_theta_K();

  SparseMatrix<double> Kl = left->getK(theta_K.segment(0, n1));
  SparseMatrix<double> Kr = right->getK(theta_K.segment(n1, n2));

  return kronecker(Kl, Kr);
}

SparseMatrix<double> Tensor_prod::get_dK(int index, const VectorXd& alpha) const {
  int n1 = left->get_n_theta_K();
  int n2 = right->get_n_theta_K();

  if (index < n1) {
    SparseMatrix<double> dk_l = left->get_dK_by_index(index);
    SparseMatrix<double> K_r = right->getK();
    return kronecker(dk_l, K_r);
  } else {
    SparseMatrix<double> dk_r = right->get_dK_by_index(index - n1);
    SparseMatrix<double> K_l = left->getK();
    return kronecker(K_l, dk_r);
  }
}

// after set parameters, update K and K_rep
void Tensor_prod::update_each_iter() {
  int n1 = left->get_n_theta_K();
  int n2 = right->get_n_theta_K();
  // setK for left and right
  left->set_theta_K(theta_K.segment(0, n1));
  right->set_theta_K(theta_K.segment(n1, n2));
// std::cout << "Theta K = " << theta_K << std::endl;

  SparseMatrix<double> K1 = left->getK();
  SparseMatrix<double> K2 = right->getK();
  K = kronecker(K1, K2);
  for (int i=0; i < n_rep; i++)
    setSparseBlock(&K_rep, i*V_size, i*W_size, K);

  std::cout << "Finish update" << std::endl;
  // SparseMatrix<double> Q = K.transpose() * K;

  // dKtx_t = dKt_t %x% Kx
  // dKtx_x = kt %x% dKx_x

  // for now, assume 1 parameter for each only
  // dK  = kronecker(left->get_dK(0), K2);
  // dK2 = kronecker(K1, right->get_dK(0));
}

