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
double time = 0;
auto timer_computeg = std::chrono::steady_clock::now();
// std::cout << "update K now" << std::endl;
  first->update_K(theta_K.segment(0, n_theta_1));
  second->update_K(theta_K.segment(n_theta_1, n_theta_2));

  // use Eigen kronecker product
  KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->getK(), second->getK());

  kroneckerEigen.evalTo(K);
time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timer_computeg).count();
// std::cout << "size and time for kronecker product is " << K.rows() << " " << K.cols() << " " << time << std::endl;
}

void Tensor_prod::update_dK(const VectorXd& theta_K) {
  VectorXd theta_K_1 = theta_K.segment(0, n_theta_1);
  VectorXd theta_K_2 = theta_K.segment(n_theta_1, n_theta_2);
  // assume K is already updated!!
  first->update_dK(theta_K_1);
  second->update_dK(theta_K_2);
  for (int index=0; index < n_theta_1 + n_theta_2; index++) {
    if (index < n_theta_1) {
      // return kroneckerEigen(dK_1, K_2);
      KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->get_dK()[index], second->getK());
      kroneckerEigen.evalTo(dK[index]);
    } else {
      // return kroneckerEigen(K_1, dK_2);
      KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->getK(), second->get_dK()[index - n_theta_1]);
      kroneckerEigen.evalTo(dK[index]);
    }
  }
}

