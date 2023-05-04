// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x

#include "../operator.h"
#include "MatrixAlgebra.h"

// tensor product for the C G class
Tensor_prod::Tensor_prod(const Rcpp::List& operator_list):
  Operator(operator_list),
  first (OperatorFactory::create(operator_list["first"])),
  second (OperatorFactory::create(operator_list["second"])),
  n_theta_1 (first->get_n_theta_K()),
  n_theta_2 (second->get_n_theta_K())
{}

SparseMatrix<double> Tensor_prod::getK(const VectorXd& theta_K) const {
  SparseMatrix<double> K1 = first->getK(theta_K.segment(0, n_theta_1));
  SparseMatrix<double> K2 = second->getK(theta_K.segment(n_theta_1, n_theta_2));

  return kroneckerEigen(K1, K2);
}

SparseMatrix<double> Tensor_prod::get_dK(int index, const VectorXd& theta_K) const {
  VectorXd theta_K_1 = theta_K.segment(0, n_theta_1);
  VectorXd theta_K_2 = theta_K.segment(n_theta_1, n_theta_2);

  if (index < n_theta_1) {
    SparseMatrix<double> dK_1 = first->get_dK(index, theta_K_1);
    SparseMatrix<double> K_2 = second->getK(theta_K_2);
    return kroneckerEigen(dK_1, K_2);
  } else {
    SparseMatrix<double> dK_2 = second->get_dK(index - n_theta_2, theta_K_2);
    SparseMatrix<double> K_1 = first->getK(theta_K_1);
    return kroneckerEigen(K_1, dK_2);
  }
}


// ------ iid model -----
Iid::Iid(const Rcpp::List& operator_list):
  Operator(operator_list),
  I(Rcpp::as< SparseMatrix<double,0,int> > (operator_list["K"])) {}

SparseMatrix<double> Iid::getK(const VectorXd& alpha) const {return I;};
SparseMatrix<double> Iid::get_dK(int index, const VectorXd& alpha) const {return 0*I;};

