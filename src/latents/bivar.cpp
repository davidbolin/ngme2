#include "../operator.h"
#include "MatrixAlgebra.h"

// ------- bivariate model -------
Bivar::Bivar(const Rcpp::List& operator_list):
  Operator(operator_list),
  first (OperatorFactory::create(operator_list["first"])),
  second (OperatorFactory::create(operator_list["second"])),
  n_theta_1 (first->get_n_theta_K()),
  n_theta_2 (second->get_n_theta_K()),
  n (h.size() / 2),
  share_param (Rcpp::as<bool> (operator_list["share_param"]))
{}

void Bivar::update_K(const VectorXd& theta_K) {
  double theta = theta_K(0);
  double rho = theta_K(1);
  VectorXd theta_K1 = theta_K.segment(2, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(2 + n_theta_1, n_theta_2);

  Matrix2d D = getD(theta, rho);
  first->update_K(theta_K1);
  second->update_K(theta_K2);

  SparseMatrix<double> K00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * second->getK();
  SparseMatrix<double> K10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * second->getK();

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);
}

void Bivar::update_dK(const VectorXd& theta_K) {
  double theta = theta_K(0);
  double rho = theta_K(1);
  VectorXd theta_K1 = theta_K.segment(2, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(2 + n_theta_1, n_theta_2);
  // assume K is updated!!!
    // first->update_K(theta_K1);
    // second->update_K(theta_K2);
  first->update_dK(theta_K1);
  second->update_dK(theta_K2);

  for (int index=0; index < n_theta_K; index++) {
    dK[index].setZero();
    if (index == 0 || index == 1) {
      // d_theta K = Dtheta * K
      SparseMatrix<double> K1 = first->getK();
      SparseMatrix<double> K2 = second->getK();
      Matrix2d dD;
      if (index == 0)
        dD = get_dD_theta(theta, rho);
      else
        dD = get_dD_rho(theta, rho);

      SparseMatrix<double> dK00 = VectorXd::Constant(n, dD(0,0)).asDiagonal() * K1;
      SparseMatrix<double> dK01 = VectorXd::Constant(n, dD(0,1)).asDiagonal() * K2;
      SparseMatrix<double> dK10 = VectorXd::Constant(n, dD(1,0)).asDiagonal() * K1;
      SparseMatrix<double> dK11 = VectorXd::Constant(n, dD(1,1)).asDiagonal() * K2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
    } else if (!share_param && index < 2 + n_theta_1) {
      Matrix2d D = getD(theta, rho);
      SparseMatrix<double> dK00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * first->get_dK()[index-2];
      SparseMatrix<double> dK10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * first->get_dK()[index-2];
      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], n, 0, dK10);
    } else if (!share_param) {
      Matrix2d D = getD(theta, rho);
      int index_in_second = index - 2 - n_theta_1;
      SparseMatrix<double> dK01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * second->get_dK()[index_in_second];
      SparseMatrix<double> dK11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * second->get_dK()[index_in_second];
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, n, dK11);
    } else {
      // share param case
      Matrix2d D = getD(theta, rho);
      SparseMatrix<double> dK00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * first->get_dK()[index-2];
      SparseMatrix<double> dK10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * first->get_dK()[index-2];
      SparseMatrix<double> dK2 = second->get_dK()[index-2];
      SparseMatrix<double> dK01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * dK2;
      SparseMatrix<double> dK11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * dK2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
    }
  }
}

Matrix2d Bivar::getD(double theta, double rho) const {
  Matrix2d D;
  D(0,0) = cos(theta) + rho*sin(theta);
  D(0,1) = -sin(theta) * pow(1+pow(rho,2),0.5);
  D(1,0) = sin(theta) - rho*cos(theta);
  D(1,1) = cos(theta) * pow(1+pow(rho,2),0.5);
  return D;
}

Matrix2d Bivar::get_dD_theta(double theta, double rho) const {
  Matrix2d Dtheta;
  Dtheta(0,0) = -sin(theta) + rho*cos(theta);
  Dtheta(0,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  Dtheta(1,0) = cos(theta) + rho*sin(theta);
  Dtheta(1,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta;
}

Matrix2d Bivar::get_dD_rho(double theta, double rho) const {
  Matrix2d Drho;
  Drho(0,0) = sin(theta);
  Drho(0,1) = -sin(theta)*rho*pow(1+pow(rho,2),-0.5);
  Drho(1,0) = -cos(theta);
  Drho(1,1) = cos(theta)*rho*pow(1+pow(rho,2),-0.5);
  return Drho;
}

Matrix2d Bivar::get_dD2_theta(double theta, double rho) const {
  Matrix2d Dtheta2;
  Dtheta2(0,0) = -cos(theta) - rho*sin(theta);
  Dtheta2(0,1) = sin(theta)*pow(1+pow(rho,2),0.5);
  Dtheta2(1,0) = -sin(theta) + rho*cos(theta);
  Dtheta2(1,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta2;
}

Matrix2d Bivar::get_dD2_rho(double theta, double rho) const {
  MatrixXd Drho2(2,2);
  Drho2(0,0) = 0;
  Drho2(0,1) = -sin(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(0,1)+= sin(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  Drho2(1,0) = 0;
  Drho2(1,1) = cos(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(1,1) -= cos(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  return Drho2;
}