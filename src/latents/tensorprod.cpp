// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x


#include "../operator.h"
#include "MatrixAlgebra.h"

// tensor product for the C G class
Tensor_prod::Tensor_prod(const Rcpp::List& operator_list):
  Operator(operator_list) {
  std::cout << "start init" << std::endl;

  Rcpp::List left_list = operator_list["left"];
  std::string model = Rcpp::as<std::string>(left_list["model"]);
  left = OperatorFactory::create(left_list);

  Rcpp::List right_list = operator_list["right"];
  model = Rcpp::as<std::string>(right_list["model"]);
  right = OperatorFactory::create(right_list);

}

SparseMatrix<double> Tensor_prod::getK(const VectorXd& theta_K) const {
  SparseMatrix<double> Kl = left->getK(theta_K.segment(0, n_theta_l));
  SparseMatrix<double> Kr = right->getK(theta_K.segment(n_theta_l, n_theta_r));

  return kronecker(Kl, Kr);
}

SparseMatrix<double> Tensor_prod::get_dK(int index, const VectorXd& theta_K) const {
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

Iid::Iid(const Rcpp::List& operator_list):
  Operator(operator_list),
  I(Rcpp::as< SparseMatrix<double,0,int> > (operator_list["K"])) {}

SparseMatrix<double> Iid::getK(const VectorXd& alpha) const {return I;};
SparseMatrix<double> Iid::get_dK(int index, const VectorXd& alpha) const {return 0*I;};

