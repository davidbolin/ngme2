// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x

#include <unsupported/Eigen/KroneckerProduct>
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

void Tensor_prod::update_K(const VectorXd& theta_K) {
  // report the time for this function
// double time = 0;
// auto timer_computeg = std::chrono::steady_clock::now();
  first->update_K(theta_K.segment(0, n_theta_1));
  second->update_K(theta_K.segment(n_theta_1, n_theta_2));

  // use Eigen kronecker product
  KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->getK(), second->getK());

  kroneckerEigen.evalTo(K);
// time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timer_computeg).count();
// std::cout << "size and time for kronecker product is " << K.rows() << " " << K.cols() << " " << time << std::endl;
}

SparseMatrix<double> Tensor_prod::get_dK(int index, const VectorXd& theta_K) const {
  VectorXd theta_K_1 = theta_K.segment(0, n_theta_1);
  VectorXd theta_K_2 = theta_K.segment(n_theta_1, n_theta_2);

  // to-do
  // if (index < n_theta_1) {
  //   SparseMatrix<double> dK_1 = first->get_dK(index, theta_K_1);
  //   SparseMatrix<double> K_2 = second->getK(theta_K_2);
  //   return kroneckerEigen(dK_1, K_2);
  // } else {
  //   SparseMatrix<double> dK_2 = second->get_dK(index - n_theta_2, theta_K_2);
  //   SparseMatrix<double> K_1 = first->getK(theta_K_1);
  //   return kroneckerEigen(K_1, dK_2);
  // }
}


// ------ iid model -----
Iid::Iid(const Rcpp::List& operator_list):
  Operator(operator_list),
  I(Rcpp::as< SparseMatrix<double,0,int> > (operator_list["K"])) {}

void Iid::update_K(const VectorXd& alpha) {}
SparseMatrix<double> Iid::get_dK(int index, const VectorXd& alpha) const {return 0*I;}

