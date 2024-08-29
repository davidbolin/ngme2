#include "../operator.h"
#include "MatrixAlgebra.h"
// #include <chrono>

// ------- bivariate model -------
Bivar::Bivar(const Rcpp::List& operator_list):
  Operator(operator_list),
  first (OperatorFactory::create(operator_list["first"])),
  second (OperatorFactory::create(operator_list["second"])),
  n_theta_1 (first->get_n_theta_K()),
  n_theta_2 (second->get_n_theta_K()),
  n (h.size() / 2),
  share_param (Rcpp::as<bool> (operator_list["share_param"])),
  fix_bv_theta (Rcpp::as<bool> (operator_list["fix_bv_theta"]))
{}

void Bivar::update_K(const VectorXd& theta_K) {
  double eta = theta_K(0);
  double theta = std::atan(eta);
  if (fix_bv_theta) { theta = 0; }

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

// assume K is updated!!!
void Bivar::update_dK(const VectorXd& theta_K) {
  double eta = theta_K(0);
  double theta = std::atan(eta);
  if (fix_bv_theta) { theta = 0; }

  double rho = theta_K(1);
  VectorXd theta_K1 = theta_K.segment(2, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(2 + n_theta_1, n_theta_2);

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
      if (index == 0) { // theta
        // d0 is dtheta / deta
      // double theta = std::atan(eta) / 2;
        double d0 = 1 / (2 * (1 + pow(eta, 2)));
        dD = d0 * get_dD_theta(theta, rho);
      }
      else { // rho
        dD = get_dD_rho(theta, rho);
      }
// std::cout << "dD = " << dD << std::endl;
      SparseMatrix<double> dK00 = VectorXd::Constant(n, dD(0,0)).asDiagonal() * K1;
      SparseMatrix<double> dK01 = VectorXd::Constant(n, dD(0,1)).asDiagonal() * K2;
      SparseMatrix<double> dK10 = VectorXd::Constant(n, dD(1,0)).asDiagonal() * K1;
      SparseMatrix<double> dK11 = VectorXd::Constant(n, dD(1,1)).asDiagonal() * K2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
// std::cout << "dK[" << index << "] = " << dK[index] << std::endl;
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


// ------- bivariate model (theta=0, Gaussian noise) -------
Bivar_normal_ope::Bivar_normal_ope(const Rcpp::List& operator_list):
  Operator(operator_list),
  first (OperatorFactory::create(operator_list["first"])),
  second (OperatorFactory::create(operator_list["second"])),
  n_theta_1 (first->get_n_theta_K()),
  n_theta_2 (second->get_n_theta_K()),
  n (h.size() / 2),
  share_param (Rcpp::as<bool> (operator_list["share_param"])),
  fix_bv_theta (Rcpp::as<bool> (operator_list["fix_bv_theta"]))
{}

void Bivar_normal_ope::update_K(const VectorXd& theta_K) {
  double rho = theta_K(0);
  double c1 = exp(theta_K(1));
  double c2 = exp(theta_K(2));
  VectorXd theta_K1 = theta_K.segment(3, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(3 + n_theta_1, n_theta_2);
  Matrix2d D = getD(0, rho);
  first->update_K(theta_K1);
  second->update_K(theta_K2);

  SparseMatrix<double> K00 = c1 * VectorXd::Constant(n, D(0,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K01 = c2 * VectorXd::Constant(n, D(0,1)).asDiagonal() * second->getK();
  SparseMatrix<double> K10 = c1 * VectorXd::Constant(n, D(1,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K11 = c2 * VectorXd::Constant(n, D(1,1)).asDiagonal() * second->getK();

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);
}

// assume K is updated!!!
void Bivar_normal_ope::update_dK(const VectorXd& theta_K) {
  double rho = theta_K(0);
  VectorXd theta_K1 = theta_K.segment(1, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(1 + n_theta_1, n_theta_2);

  // first->update_K(theta_K1);
  // second->update_K(theta_K2);
  first->update_dK(theta_K1);
  second->update_dK(theta_K2);

  for (int index=0; index < n_theta_K; index++) {
    dK[index].setZero();
    if (index == 0) {
      // d_theta K = Dtheta * K
      SparseMatrix<double> K1 = first->getK();
      SparseMatrix<double> K2 = second->getK();
      Matrix2d dD;
      dD = get_dD_rho(0, rho);

      SparseMatrix<double> dK00 = VectorXd::Constant(n, dD(0,0)).asDiagonal() * K1;
      SparseMatrix<double> dK01 = VectorXd::Constant(n, dD(0,1)).asDiagonal() * K2;
      SparseMatrix<double> dK10 = VectorXd::Constant(n, dD(1,0)).asDiagonal() * K1;
      SparseMatrix<double> dK11 = VectorXd::Constant(n, dD(1,1)).asDiagonal() * K2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);

    } else if (!share_param && index < 1 + n_theta_1) {
      Matrix2d D = getD(0, rho);
      SparseMatrix<double> dK00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * first->get_dK()[index-1];
      SparseMatrix<double> dK10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * first->get_dK()[index-1];
      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], n, 0, dK10);
    } else if (!share_param) {
      Matrix2d D = getD(0, rho);
      int index_in_second = index - 1 - n_theta_1;
      SparseMatrix<double> dK01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * second->get_dK()[index_in_second];
      SparseMatrix<double> dK11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * second->get_dK()[index_in_second];
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, n, dK11);
    } else {
      // share param case
      Matrix2d D = getD(0, rho);
      SparseMatrix<double> dK00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * first->get_dK()[index-1];
      SparseMatrix<double> dK10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * first->get_dK()[index-1];
      SparseMatrix<double> dK2 = second->get_dK()[index-1];
      SparseMatrix<double> dK01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * dK2;
      SparseMatrix<double> dK11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * dK2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
    }
  }
}

Matrix2d Bivar_normal_ope::getD(double theta, double rho) const {
  Matrix2d D;
  D(0,0) = cos(theta) + rho*sin(theta);
  D(0,1) = -sin(theta) * pow(1+pow(rho,2),0.5);
  D(1,0) = sin(theta) - rho*cos(theta);
  D(1,1) = cos(theta) * pow(1+pow(rho,2),0.5);
  return D;
}

Matrix2d Bivar_normal_ope::get_dD_theta(double theta, double rho) const {
  Matrix2d Dtheta;
  Dtheta(0,0) = -sin(theta) + rho*cos(theta);
  Dtheta(0,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  Dtheta(1,0) = cos(theta) + rho*sin(theta);
  Dtheta(1,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta;
}

Matrix2d Bivar_normal_ope::get_dD_rho(double theta, double rho) const {
  Matrix2d Drho;
  Drho(0,0) = sin(theta);
  Drho(0,1) = -sin(theta)*rho*pow(1+pow(rho,2),-0.5);
  Drho(1,0) = -cos(theta);
  Drho(1,1) = cos(theta)*rho*pow(1+pow(rho,2),-0.5);
  return Drho;
}

Matrix2d Bivar_normal_ope::get_dD2_theta(double theta, double rho) const {
  Matrix2d Dtheta2;
  Dtheta2(0,0) = -cos(theta) - rho*sin(theta);
  Dtheta2(0,1) = sin(theta)*pow(1+pow(rho,2),0.5);
  Dtheta2(1,0) = -sin(theta) + rho*cos(theta);
  Dtheta2(1,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta2;
}

Matrix2d Bivar_normal_ope::get_dD2_rho(double theta, double rho) const {
  MatrixXd Drho2(2,2);
  Drho2(0,0) = 0;
  Drho2(0,1) = -sin(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(0,1)+= sin(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  Drho2(1,0) = 0;
  Drho2(1,1) = cos(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(1,1) -= cos(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  return Drho2;
}



// ------- bivariate model only for matern normal model
bv_matern_normal::bv_matern_normal(const Rcpp::List& operator_list):
  Operator(operator_list),
  first (std::make_shared<Matern>(operator_list["first"])),
  second (std::make_shared<Matern>(operator_list["second"])),
  n_theta_1 (first->get_n_theta_K()),
  n_theta_2 (second->get_n_theta_K()),
  n (h.size() / 2),
  share_param (Rcpp::as<bool> (operator_list["share_param"])),
  fix_bv_theta (Rcpp::as<bool> (operator_list["fix_bv_theta"])),
  dim (Rcpp::as<int> (operator_list["dim"])),
  alpha1 (first->get_alpha()),
  alpha2 (second->get_alpha()),
  nu1 (alpha1 - 0.5*dim),
  nu2 (alpha2 - 0.5*dim)
{}

void bv_matern_normal::update_K(const VectorXd& theta_K) {
  double rho = theta_K(0);
  double sd1 = exp(theta_K(1));
  double sd2 = exp(theta_K(2));
  VectorXd theta_K1 = theta_K.segment(3, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(3 + n_theta_1, n_theta_2);
  Matrix2d D = getD(0, rho);
  first->update_K(theta_K1);
  second->update_K(theta_K2);

  double kappa1 = exp(theta_K1[0]);
  double kappa2 = exp(theta_K2[0]);
// c1 = sqrt(gamma(nu1) / (gamma(alpha1))) / (kappa1^(nu1) * (4*pi)^(d/4) * sd1)

  // compute c1, and c2
  double c1 = sqrt(tgamma(nu1) / tgamma(alpha1)) / pow(kappa1, nu1) / pow(4*M_PI, 0.25 * dim) / sd1;
  double c2 = sqrt(tgamma(nu2) / tgamma(alpha2)) / pow(kappa2, nu2) / pow(4*M_PI, 0.25 * dim) / sd2;

  SparseMatrix<double> K00 = c1 * VectorXd::Constant(n, D(0,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K01 = c2 * VectorXd::Constant(n, D(0,1)).asDiagonal() * second->getK();
  SparseMatrix<double> K10 = c1 * VectorXd::Constant(n, D(1,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K11 = c2 * VectorXd::Constant(n, D(1,1)).asDiagonal() * second->getK();

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);
}

void bv_matern_normal::update_dK(const VectorXd& theta_K) {
  // no supported
}

Matrix2d bv_matern_normal::getD(double theta, double rho) const {
  Matrix2d D;
  D(0,0) = cos(theta) + rho*sin(theta);
  D(0,1) = -sin(theta) * pow(1+pow(rho,2),0.5);
  D(1,0) = sin(theta) - rho*cos(theta);
  D(1,1) = cos(theta) * pow(1+pow(rho,2),0.5);
  return D;
}

Matrix2d bv_matern_normal::get_dD_theta(double theta, double rho) const {
  Matrix2d Dtheta;
  Dtheta(0,0) = -sin(theta) + rho*cos(theta);
  Dtheta(0,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  Dtheta(1,0) = cos(theta) + rho*sin(theta);
  Dtheta(1,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta;
}

Matrix2d bv_matern_normal::get_dD_rho(double theta, double rho) const {
  Matrix2d Drho;
  Drho(0,0) = sin(theta);
  Drho(0,1) = -sin(theta)*rho*pow(1+pow(rho,2),-0.5);
  Drho(1,0) = -cos(theta);
  Drho(1,1) = cos(theta)*rho*pow(1+pow(rho,2),-0.5);
  return Drho;
}

Matrix2d bv_matern_normal::get_dD2_theta(double theta, double rho) const {
  Matrix2d Dtheta2;
  Dtheta2(0,0) = -cos(theta) - rho*sin(theta);
  Dtheta2(0,1) = sin(theta)*pow(1+pow(rho,2),0.5);
  Dtheta2(1,0) = -sin(theta) + rho*cos(theta);
  Dtheta2(1,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta2;
}

Matrix2d bv_matern_normal::get_dD2_rho(double theta, double rho) const {
  MatrixXd Drho2(2,2);
  Drho2(0,0) = 0;
  Drho2(0,1) = -sin(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(0,1)+= sin(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  Drho2(1,0) = 0;
  Drho2(1,1) = cos(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(1,1) -= cos(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  return Drho2;
}



// ------- bivariate model only for matern normal model
bv_matern_nig::bv_matern_nig(const Rcpp::List& operator_list):
  Operator(operator_list),
  first (std::make_shared<Matern>(operator_list["first"])),
  second (std::make_shared<Matern>(operator_list["second"])),
  n_theta_1 (first->get_n_theta_K()),
  n_theta_2 (second->get_n_theta_K()),
  n (h.size() / 2),
  share_param (Rcpp::as<bool> (operator_list["share_param"])),
  fix_bv_theta (Rcpp::as<bool> (operator_list["fix_bv_theta"])),
  dim (Rcpp::as<int> (operator_list["dim"])),
  alpha1 (first->get_alpha()),
  alpha2 (second->get_alpha()),
  nu1 (alpha1 - 0.5*dim),
  nu2 (alpha2 - 0.5*dim)
{}

void bv_matern_nig::update_K(const VectorXd& theta_K) {
  double theta = theta_K(0);
  double rho = theta_K(1);
  double sd1 = exp(theta_K(2));
  double sd2 = exp(theta_K(3));
  VectorXd theta_K1 = theta_K.segment(4, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(4 + n_theta_1, n_theta_2);
  Matrix2d D = getD(theta, rho);
  first->update_K(theta_K1);
  second->update_K(theta_K2);

  double kappa1 = exp(theta_K1[0]);
  double kappa2 = exp(theta_K2[0]);
// c1 = sqrt(gamma(nu1) / (gamma(alpha1))) / (kappa1^(nu1) * (4*pi)^(d/4) * sd1)

  // compute c1, and c2
  double c1 = sqrt(tgamma(nu1) / tgamma(alpha1)) / pow(kappa1, nu1) / pow(4*M_PI, 0.25 * dim) / sd1;
  double c2 = sqrt(tgamma(nu2) / tgamma(alpha2)) / pow(kappa2, nu2) / pow(4*M_PI, 0.25 * dim) / sd2;

  SparseMatrix<double> K00 = c1 * VectorXd::Constant(n, D(0,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K01 = c2 * VectorXd::Constant(n, D(0,1)).asDiagonal() * second->getK();
  SparseMatrix<double> K10 = c1 * VectorXd::Constant(n, D(1,0)).asDiagonal() * first->getK();
  SparseMatrix<double> K11 = c2 * VectorXd::Constant(n, D(1,1)).asDiagonal() * second->getK();

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);
}

void bv_matern_nig::update_dK(const VectorXd& theta_K) {
  // no supported
}

Matrix2d bv_matern_nig::getD(double theta, double rho) const {
  Matrix2d D;
  D(0,0) = cos(theta) + rho*sin(theta);
  D(0,1) = -sin(theta) * pow(1+pow(rho,2),0.5);
  D(1,0) = sin(theta) - rho*cos(theta);
  D(1,1) = cos(theta) * pow(1+pow(rho,2),0.5);
  return D;
}

Matrix2d bv_matern_nig::get_dD_theta(double theta, double rho) const {
  Matrix2d Dtheta;
  Dtheta(0,0) = -sin(theta) + rho*cos(theta);
  Dtheta(0,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  Dtheta(1,0) = cos(theta) + rho*sin(theta);
  Dtheta(1,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta;
}

Matrix2d bv_matern_nig::get_dD_rho(double theta, double rho) const {
  Matrix2d Drho;
  Drho(0,0) = sin(theta);
  Drho(0,1) = -sin(theta)*rho*pow(1+pow(rho,2),-0.5);
  Drho(1,0) = -cos(theta);
  Drho(1,1) = cos(theta)*rho*pow(1+pow(rho,2),-0.5);
  return Drho;
}

Matrix2d bv_matern_nig::get_dD2_theta(double theta, double rho) const {
  Matrix2d Dtheta2;
  Dtheta2(0,0) = -cos(theta) - rho*sin(theta);
  Dtheta2(0,1) = sin(theta)*pow(1+pow(rho,2),0.5);
  Dtheta2(1,0) = -sin(theta) + rho*cos(theta);
  Dtheta2(1,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  return Dtheta2;
}

Matrix2d bv_matern_nig::get_dD2_rho(double theta, double rho) const {
  MatrixXd Drho2(2,2);
  Drho2(0,0) = 0;
  Drho2(0,1) = -sin(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(0,1)+= sin(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  Drho2(1,0) = 0;
  Drho2(1,1) = cos(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(1,1) -= cos(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  return Drho2;
}



