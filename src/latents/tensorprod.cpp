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

void Tensor_prod::build_KZ(const VectorXd& theta_K) {
  // report the time for this function
double time = 0;
auto timer_computeg = std::chrono::steady_clock::now();
// std::cout << "update K now" << std::endl;
  first->build_KZ(theta_K.segment(0, n_theta_1));
// std::cout << "update K now1" << std::endl;
  second->build_KZ(theta_K.segment(n_theta_1, n_theta_2));
// std::cout << "new K size = " << first->getK().rows() * second->getK().rows() << " " << first->getK().cols() * second->getK().cols() << std::endl;

  // use Eigen kronecker product
  KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->getK(), second->getK());

  kroneckerEigen.evalTo(K);
// std::cout << "update K now3" << std::endl;

time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timer_computeg).count();
// std::cout << "size and time for kronecker product is " << K.rows() << " " << K.cols() << " " << time << std::endl;
  // Z remains identity (set once in Operator base)
}

bool Tensor_prod::update_dKdZ(const VectorXd& theta_K) {
  VectorXd theta_K_1 = theta_K.segment(0, n_theta_1);
  VectorXd theta_K_2 = theta_K.segment(n_theta_1, n_theta_2);
  // assume K is already updated!!
  first->update_dKdZ(theta_K_1);
  second->update_dKdZ(theta_K_2);
  for (int index=0; index < n_theta_1 + n_theta_2; index++) {
    if (index < n_theta_1) {
      // return kroneckerEigen(dK_1, K_2);
      KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->get_dK(index), second->getK());
      kroneckerEigen.evalTo(dK[index]);
    } else {
      // return kroneckerEigen(K_1, dK_2);
      KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(first->getK(), second->get_dK(index - n_theta_1));
      kroneckerEigen.evalTo(dK[index]);
    }
  }
  return true;
}



// Non-separable Space-time model

Spacetime::Spacetime(const Rcpp::List& operator_list):
  Operator(operator_list),
  Ct_diag (Rcpp::as<VectorXd> (operator_list["Ct_diag"])),
  Cs_diag (Rcpp::as<VectorXd> (operator_list["Cs_diag"])),
  BtCs (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["BtCs"])),
  Gs (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Gs"])),
  Ct (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Ct"])),
  Cs (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Cs"])),
  Bx (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Bx"])),
  By (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["By"])),
  S (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["S"])),
  Bs (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Bs"])),
  Hxx (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Hxx"])),
  Hyy (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Hyy"])),
  Hxy (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Hxy"])),
  Hyx (Rcpp::as<SparseMatrix<double, 0, int>> (operator_list["Hyx"])),
  // B_gamma_x (Rcpp::as<MatrixXd> (operator_list["B_gamma_x"])),
  // B_gamma_y (Rcpp::as<MatrixXd> (operator_list["B_gamma_y"])),
  B_gamma_x_list_input (Rcpp::as<Rcpp::List> (operator_list["B_gamma_x_list"])),
  B_gamma_y_list_input (Rcpp::as<Rcpp::List> (operator_list["B_gamma_y_list"])),
  theta_gamma_x (Rcpp::as<VectorXd> (operator_list["theta_gamma_x"])),
  theta_gamma_y (Rcpp::as<VectorXd> (operator_list["theta_gamma_y"])),
  n_theta_gamma_x (Rcpp::as<int> (operator_list["n_theta_gamma_x"])),
  n_theta_gamma_y (Rcpp::as<int> (operator_list["n_theta_gamma_y"])),
  lambda (Rcpp::as<double> (operator_list["lambda"])),
  alpha (Rcpp::as<double> (operator_list["alpha"])),
  method (Rcpp::as<string> (operator_list["method"])),
  stabilization (Rcpp::as<bool> (operator_list["stabilization"])),
  fix_gamma (Rcpp::as<bool> (operator_list["fix_gamma"])),
  shared_theta_gamma (Rcpp::as<bool> (operator_list["shared_theta_gamma"])),
  nt (Rcpp::as<int> (operator_list["nt"])),
  B_gamma_x_list (nt-1),
  B_gamma_y_list (nt-1)
{
  // turn B_gamma_x_list_input into B_gamma_x_list
  if (!fix_gamma) {
    for (int i = 0; i < nt-1; i++) {
      B_gamma_x_list[i] = Rcpp::as<MatrixXd> (B_gamma_x_list_input[i]);
      B_gamma_y_list[i] = Rcpp::as<MatrixXd> (B_gamma_y_list_input[i]);
    }
  }
}

void Spacetime::build_KZ(const VectorXd& theta_K) {
  K.setZero();
  double c = exp(theta_K[0]);
  double kappa = exp(theta_K[1]);

  if (!fix_gamma) {
    if (shared_theta_gamma) {
      theta_gamma_x = theta_K.segment(2, n_theta_gamma_x);
      theta_gamma_y = theta_gamma_x;  // Use same parameter for both x and y
    } else {
      theta_gamma_x = theta_K.segment(2, n_theta_gamma_x);
      theta_gamma_y = theta_K.segment(2 + n_theta_gamma_x, n_theta_gamma_y);
    }
  
    std::vector<VectorXd> gamma_x_list (nt-1);
    std::vector<VectorXd> gamma_y_list (nt-1);
    for (int i = 0; i < nt-1; i++) {
      gamma_x_list[i] = B_gamma_x_list[i] * theta_gamma_x;
      gamma_y_list[i] = B_gamma_y_list[i] * theta_gamma_y;
    }

    // Build Bs_list, S_list, Ls_list
    std::vector<SparseMatrix<double, 0, int>> Bs_list (nt-1);
    std::vector<SparseMatrix<double, 0, int>> S_list (nt-1);
    std::vector<SparseMatrix<double, 0, int>> Ls_list (nt-1);
    for (int i = 0; i < nt-1; i++) {
      Bs_list[i] = gamma_x_list[i].asDiagonal() * Bx * gamma_x_list[i].asDiagonal() + 
        gamma_y_list[i].asDiagonal() * By * gamma_y_list[i].asDiagonal();

      Ls_list[i] = kappa*kappa * Cs + lambda * Gs + Bs_list[i];

      if (alpha == 4) {
        Ls_list[i] = Ls_list[i] * Cs_diag.cwiseInverse().asDiagonal() * Ls_list[i].transpose();
      }

      if (stabilization) {
        // Compute S_list[i]
        if (gamma_x_list[i].norm() < 1e-8 && gamma_y_list[i].norm() < 1e-8) {
          S_list[i] = SparseMatrix<double, 0, int>(Ls_list[i].rows(), Ls_list[i].cols());
          S_list[i].setZero();
        } else {
          VectorXd gamma_xx = gamma_x_list[i].array().square();
          VectorXd gamma_yy = gamma_y_list[i].array().square();
          VectorXd gamma_xy = gamma_x_list[i].array() * gamma_y_list[i].array();

          S_list[i] = gamma_xx.asDiagonal() * Hxx * gamma_xx.asDiagonal() + 
            gamma_yy.asDiagonal() * Hyy * gamma_yy.asDiagonal() + 
            gamma_xy.asDiagonal() * (Hxy + Hyx) * gamma_xy.asDiagonal();

          double gamma_norm = sqrt((gamma_x_list[i].array().square() + gamma_y_list[i].array().square()).sum());
          S_list[i] = Cs_diag.asDiagonal() * S_list[i] / gamma_norm;
        }
        // Add S_list[i] to Ls_list[i]
        Ls_list[i] = Ls_list[i] + S_list[i];
      }
    }

    for (int i = 1; i < Ct.rows(); ++i) {
      setSparseBlock(&K, i * Ls_list[i-1].rows(), i * Ls_list[i-1].cols(), Ls_list[i-1]);
    }
  } else { // fix_gamma = TRUE
    SparseMatrix<double> Ls = (kappa*kappa * Cs + lambda * Gs + Bs);
    // alpha=4, L = L %*% solve(Ct %x% Cs, L) 
    if (alpha == 4) 
      Ls = Ls * Cs_diag.cwiseInverse().asDiagonal() * Ls.transpose();
    for (int i = 1; i < Ct.rows(); ++i) {
      setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls);
    }
  }
  K = BtCs + K / c;

  K = sqrt(c) * K;

  // SparseMatrix<double> Ls = (kappa*kappa * Cs + lambda * Gs + Bs);
  
  // alpha=4, L = L %*% solve(Ct %x% Cs, L) 
  // if (alpha == 4) 
  //   Ls = Ls * Cs_diag.cwiseInverse().asDiagonal() * Ls.transpose();
  // Ct is diagonal

  // if (method == "galerkin") {
  //   KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(Ct, Ls);

  //   // update K
  //   K = kroneckerEigen.eval();
  // } else if (method == "euler") {
  //   // if (stabilization) Ls = Ls + S;
  //   // Build K = bdiag(0, Ls, ..., Ls) 
  //   // 1st approach: Using Kronecker product
  //   // create diag(0, 1, 1, ..., 1) (nt-1) 1
  //   // Eigen::SparseMatrix<double> I_sparse(Ct.rows(), Ct.rows());
  //   // for (int i = 1; i < Ct.rows(); ++i) {
  //   //   I_sparse.insert(i, i) = 1;
  //   // }

  //   // KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> > kroneckerEigen(I_sparse, Ls);
  //   // K = kroneckerEigen.eval();

  //   // 2nd approach: build Bdiag(0, Ls, ..., Ls) directly
  //   K.setZero();
  //   for (int i = 0; i < nt-1; ++i) {
  //     setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls_list[i]);
  //   }
  // }
  
}

void Spacetime::update_dK(const VectorXd& theta_K) {
}
