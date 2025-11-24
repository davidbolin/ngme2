#include "../operator.h"
#include "MatrixAlgebra.h"
// #include <chrono>

// ------- bivariate model -------
Bivar::Bivar(const Rcpp::List &operator_list)
    : Operator(operator_list),
      first(OperatorFactory::create(
          Rcpp::as<Rcpp::List>(operator_list["first"]))),
      second(OperatorFactory::create(
          Rcpp::as<Rcpp::List>(operator_list["second"]))),
      n_theta_1(first->get_n_theta_K()), n_theta_2(second->get_n_theta_K()),
      n(h.size() / 2),
      share_param(Rcpp::as<bool>(operator_list["share_param"])),
      fix_theta(Rcpp::as<bool>(operator_list["fix_theta"])),
      use_c_param(Rcpp::as<bool>(operator_list["use_c_param"])),
      bv_theta(Rcpp::as<double>(operator_list["bv_theta"])) {}

void Bivar::build_KZ(const VectorXd &theta_K) {
  double theta;
  double rho;
  double c1 = 1.0, c2 = 1.0;
  VectorXd theta_K1, theta_K2;

  int idx = 0;
  if (fix_theta) {
    theta = bv_theta;
  } else {
    double eta = theta_K(idx++);
    theta = std::atan(eta);
  }
  rho = theta_K(idx++);

  if (use_c_param) {
    c1 = exp(theta_K(idx++));
    c2 = exp(theta_K(idx++));
  }

  theta_K1 = theta_K.segment(idx, n_theta_1);
  idx += n_theta_1;
  theta_K2 = theta_K.segment(idx, n_theta_2);

  Matrix2d D = getD(theta, rho);
  first->build_KZ(theta_K1);
  second->build_KZ(theta_K2);

  SparseMatrix<double> K00 =
      c1 * VectorXd::Constant(n, D(0, 0)).asDiagonal() * first->getK();
  SparseMatrix<double> K01 =
      c2 * VectorXd::Constant(n, D(0, 1)).asDiagonal() * second->getK();
  SparseMatrix<double> K10 =
      c1 * VectorXd::Constant(n, D(1, 0)).asDiagonal() * first->getK();
  SparseMatrix<double> K11 =
      c2 * VectorXd::Constant(n, D(1, 1)).asDiagonal() * second->getK();

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);
}

// assume K is updated!!!
bool Bivar::update_dKdZ(const VectorXd &theta_K) {
  double theta;
  double rho;
  double c1 = 1.0, c2 = 1.0;
  VectorXd theta_K1, theta_K2;

  int idx = 0;
  if (fix_theta) {
    theta = bv_theta;
  } else {
    double eta = theta_K(idx++);
    theta = std::atan(eta);
  }
  rho = theta_K(idx++);

  if (use_c_param) {
    c1 = exp(theta_K(idx++));
    c2 = exp(theta_K(idx++));
  }

  theta_K1 = theta_K.segment(idx, n_theta_1);
  idx += n_theta_1;
  theta_K2 = theta_K.segment(idx, n_theta_2);

  first->update_dKdZ(theta_K1);
  second->update_dKdZ(theta_K2);

  int offset = (fix_theta ? 0 : 1) + 1 + (use_c_param ? 2 : 0);

  for (int index = 0; index < n_theta_K; index++) {
    dK[index].setZero();

    if (!fix_theta && index == 0) {
      // d_theta K = Dtheta * K
      double eta = theta_K(0);
      double d0 = 1.0 / (1.0 + pow(eta, 2));
      Matrix2d dD = d0 * get_dD_theta(theta, rho);

      SparseMatrix<double> K1 = first->getK();
      SparseMatrix<double> K2 = second->getK();

      SparseMatrix<double> dK00 =
          c1 * VectorXd::Constant(n, dD(0, 0)).asDiagonal() * K1;
      SparseMatrix<double> dK01 =
          c2 * VectorXd::Constant(n, dD(0, 1)).asDiagonal() * K2;
      SparseMatrix<double> dK10 =
          c1 * VectorXd::Constant(n, dD(1, 0)).asDiagonal() * K1;
      SparseMatrix<double> dK11 =
          c2 * VectorXd::Constant(n, dD(1, 1)).asDiagonal() * K2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
      continue;
    }

    int rho_index = fix_theta ? 0 : 1;
    if (index == rho_index) {
      // d_rho K
      Matrix2d dD = get_dD_rho(theta, rho);

      SparseMatrix<double> K1 = first->getK();
      SparseMatrix<double> K2 = second->getK();

      SparseMatrix<double> dK00 =
          c1 * VectorXd::Constant(n, dD(0, 0)).asDiagonal() * K1;
      SparseMatrix<double> dK01 =
          c2 * VectorXd::Constant(n, dD(0, 1)).asDiagonal() * K2;
      SparseMatrix<double> dK10 =
          c1 * VectorXd::Constant(n, dD(1, 0)).asDiagonal() * K1;
      SparseMatrix<double> dK11 =
          c2 * VectorXd::Constant(n, dD(1, 1)).asDiagonal() * K2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
      continue;
    }

    if (use_c_param) {
      int c1_index = rho_index + 1;
      int c2_index = rho_index + 2;
      if (index == c1_index) {
        // d_c1 K = K (blocks 00 and 10)
        // K00 = c1 * ... * K1
        // dK00/d(log c1) = K00
        Matrix2d D = getD(theta, rho);
        SparseMatrix<double> K00 =
            c1 * VectorXd::Constant(n, D(0, 0)).asDiagonal() * first->getK();
        SparseMatrix<double> K10 =
            c1 * VectorXd::Constant(n, D(1, 0)).asDiagonal() * first->getK();
        setSparseBlock(&dK[index], 0, 0, K00);
        setSparseBlock(&dK[index], n, 0, K10);
        continue;
      }
      if (index == c2_index) {
        // d_c2 K = K (blocks 01 and 11)
        Matrix2d D = getD(theta, rho);
        SparseMatrix<double> K01 =
            c2 * VectorXd::Constant(n, D(0, 1)).asDiagonal() * second->getK();
        SparseMatrix<double> K11 =
            c2 * VectorXd::Constant(n, D(1, 1)).asDiagonal() * second->getK();
        setSparseBlock(&dK[index], 0, n, K01);
        setSparseBlock(&dK[index], n, n, K11);
        continue;
      }
    }

    // Sub-model parameters
    int sub_index = index - offset;
    if (!share_param && sub_index < n_theta_1) {
      Matrix2d D = getD(theta, rho);
      SparseMatrix<double> dK00 = c1 *
                                  VectorXd::Constant(n, D(0, 0)).asDiagonal() *
                                  first->get_dK(sub_index);
      SparseMatrix<double> dK10 = c1 *
                                  VectorXd::Constant(n, D(1, 0)).asDiagonal() *
                                  first->get_dK(sub_index);
      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], n, 0, dK10);
    } else if (!share_param) {
      Matrix2d D = getD(theta, rho);
      int index_in_second = sub_index - n_theta_1;
      SparseMatrix<double> dK01 = c2 *
                                  VectorXd::Constant(n, D(0, 1)).asDiagonal() *
                                  second->get_dK(index_in_second);
      SparseMatrix<double> dK11 = c2 *
                                  VectorXd::Constant(n, D(1, 1)).asDiagonal() *
                                  second->get_dK(index_in_second);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, n, dK11);
    } else {
      // share param case
      Matrix2d D = getD(theta, rho);
      SparseMatrix<double> dK00 = c1 *
                                  VectorXd::Constant(n, D(0, 0)).asDiagonal() *
                                  first->get_dK(sub_index);
      SparseMatrix<double> dK10 = c1 *
                                  VectorXd::Constant(n, D(1, 0)).asDiagonal() *
                                  first->get_dK(sub_index);
      SparseMatrix<double> dK2 = second->get_dK(sub_index);
      SparseMatrix<double> dK01 =
          c2 * VectorXd::Constant(n, D(0, 1)).asDiagonal() * dK2;
      SparseMatrix<double> dK11 =
          c2 * VectorXd::Constant(n, D(1, 1)).asDiagonal() * dK2;

      setSparseBlock(&dK[index], 0, 0, dK00);
      setSparseBlock(&dK[index], 0, n, dK01);
      setSparseBlock(&dK[index], n, 0, dK10);
      setSparseBlock(&dK[index], n, n, dK11);
    }
  }
  return true;
}

Matrix2d Bivar::getD(double theta, double rho) const {
  Matrix2d D;
  D(0, 0) = cos(theta) + rho * sin(theta);
  D(0, 1) = -sin(theta) * pow(1 + pow(rho, 2), 0.5);
  D(1, 0) = sin(theta) - rho * cos(theta);
  D(1, 1) = cos(theta) * pow(1 + pow(rho, 2), 0.5);
  return D;
}

Matrix2d Bivar::get_dD_theta(double theta, double rho) const {
  Matrix2d Dtheta;
  Dtheta(0, 0) = -sin(theta) + rho * cos(theta);
  Dtheta(0, 1) = -cos(theta) * pow(1 + pow(rho, 2), 0.5);
  Dtheta(1, 0) = cos(theta) + rho * sin(theta);
  Dtheta(1, 1) = -sin(theta) * pow(1 + pow(rho, 2), 0.5);
  return Dtheta;
}

Matrix2d Bivar::get_dD_rho(double theta, double rho) const {
  Matrix2d Drho;
  Drho(0, 0) = sin(theta);
  Drho(0, 1) = -sin(theta) * rho * pow(1 + pow(rho, 2), -0.5);
  Drho(1, 0) = -cos(theta);
  Drho(1, 1) = cos(theta) * rho * pow(1 + pow(rho, 2), -0.5);
  return Drho;
}

Matrix2d Bivar::get_dD2_theta(double theta, double rho) const {
  Matrix2d Dtheta2;
  Dtheta2(0, 0) = -cos(theta) - rho * sin(theta);
  Dtheta2(0, 1) = sin(theta) * pow(1 + pow(rho, 2), 0.5);
  Dtheta2(1, 0) = -sin(theta) + rho * cos(theta);
  Dtheta2(1, 1) = -cos(theta) * pow(1 + pow(rho, 2), 0.5);
  return Dtheta2;
}

Matrix2d Bivar::get_dD2_rho(double theta, double rho) const {
  MatrixXd Drho2(2, 2);
  Drho2(0, 0) = 0;
  Drho2(0, 1) = -sin(theta) * pow(1 + pow(rho, 2), -0.5);
  Drho2(0, 1) += sin(theta) * pow(rho, 2) * pow(1 + pow(rho, 2), -1.5);
  Drho2(1, 0) = 0;
  Drho2(1, 1) = cos(theta) * pow(1 + pow(rho, 2), -0.5);
  Drho2(1, 1) -= cos(theta) * pow(rho, 2) * pow(1 + pow(rho, 2), -1.5);
  return Drho2;
}

// ------- bivariate model only for matern normal model
bv_matern::bv_matern(const Rcpp::List &operator_list)
    : Operator(operator_list),
      first(std::make_shared<Matern>(
          Rcpp::as<Rcpp::List>(operator_list["first"]))),
      second(std::make_shared<Matern>(
          Rcpp::as<Rcpp::List>(operator_list["second"]))),
      n_theta_1(first->get_n_theta_K()), n_theta_2(second->get_n_theta_K()),
      n(h.size() / 2),
      share_param(Rcpp::as<bool>(operator_list["share_param"])),
      fix_theta(Rcpp::as<bool>(operator_list["fix_theta"])),
      dim(Rcpp::as<int>(operator_list["dim"])), alpha1(first->get_alpha()),
      alpha2(second->get_alpha()), nu1(alpha1 - 0.5 * dim),
      nu2(alpha2 - 0.5 * dim),
      bv_theta(Rcpp::as<double>(operator_list["bv_theta"])) {}

void bv_matern::build_KZ(const VectorXd &theta_K) {
  double theta;
  if (!fix_theta) {
    theta = atan(theta_K(0));
  } else {
    theta = bv_theta;
  }
  double rho = theta_K(1 - fix_theta);
  double sd1 = exp(theta_K(2 - fix_theta));
  double sd2 = exp(theta_K(3 - fix_theta));
  VectorXd theta_K1 = theta_K.segment(4 - fix_theta, n_theta_1);
  VectorXd theta_K2 = theta_K.segment(4 - fix_theta + n_theta_1, n_theta_2);
  Matrix2d D = getD(theta, rho);
  first->build_KZ(theta_K1);
  second->build_KZ(theta_K2);

  double kappa1 = exp(theta_K1[0]);
  double kappa2 = exp(theta_K2[0]);
  // c1 = sqrt(gamma(nu1) / (gamma(alpha1))) / (kappa1^(nu1) * (4*pi)^(d/4) *
  // sd1)

  // compute c1, and c2
  double c1 = sqrt(tgamma(nu1) / tgamma(alpha1)) / pow(kappa1, nu1) /
              pow(4 * M_PI, 0.25 * dim) / sd1;
  double c2 = sqrt(tgamma(nu2) / tgamma(alpha2)) / pow(kappa2, nu2) /
              pow(4 * M_PI, 0.25 * dim) / sd2;

  SparseMatrix<double> K00 =
      c1 * VectorXd::Constant(n, D(0, 0)).asDiagonal() * first->getK();
  SparseMatrix<double> K01 =
      c2 * VectorXd::Constant(n, D(0, 1)).asDiagonal() * second->getK();
  SparseMatrix<double> K10 =
      c1 * VectorXd::Constant(n, D(1, 0)).asDiagonal() * first->getK();
  SparseMatrix<double> K11 =
      c2 * VectorXd::Constant(n, D(1, 1)).asDiagonal() * second->getK();

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);
}

Matrix2d bv_matern::getD(double theta, double rho) const {
  Matrix2d D;
  D(0, 0) = cos(theta) + rho * sin(theta);
  D(0, 1) = -sin(theta) * pow(1 + pow(rho, 2), 0.5);
  D(1, 0) = sin(theta) - rho * cos(theta);
  D(1, 1) = cos(theta) * pow(1 + pow(rho, 2), 0.5);
  return D;
}

Matrix2d bv_matern::get_dD_theta(double theta, double rho) const {
  Matrix2d Dtheta;
  Dtheta(0, 0) = -sin(theta) + rho * cos(theta);
  Dtheta(0, 1) = -cos(theta) * pow(1 + pow(rho, 2), 0.5);
  Dtheta(1, 0) = cos(theta) + rho * sin(theta);
  Dtheta(1, 1) = -sin(theta) * pow(1 + pow(rho, 2), 0.5);
  return Dtheta;
}

Matrix2d bv_matern::get_dD_rho(double theta, double rho) const {
  Matrix2d Drho;
  Drho(0, 0) = sin(theta);
  Drho(0, 1) = -sin(theta) * rho * pow(1 + pow(rho, 2), -0.5);
  Drho(1, 0) = -cos(theta);
  Drho(1, 1) = cos(theta) * rho * pow(1 + pow(rho, 2), -0.5);
  return Drho;
}

Matrix2d bv_matern::get_dD2_theta(double theta, double rho) const {
  Matrix2d Dtheta2;
  Dtheta2(0, 0) = -cos(theta) - rho * sin(theta);
  Dtheta2(0, 1) = sin(theta) * pow(1 + pow(rho, 2), 0.5);
  Dtheta2(1, 0) = -sin(theta) + rho * cos(theta);
  Dtheta2(1, 1) = -cos(theta) * pow(1 + pow(rho, 2), 0.5);
  return Dtheta2;
}

Matrix2d bv_matern::get_dD2_rho(double theta, double rho) const {
  MatrixXd Drho2(2, 2);
  Drho2(0, 0) = 0;
  Drho2(0, 1) = -sin(theta) * pow(1 + pow(rho, 2), -0.5);
  Drho2(0, 1) += sin(theta) * pow(rho, 2) * pow(1 + pow(rho, 2), -1.5);
  Drho2(1, 0) = 0;
  Drho2(1, 1) = cos(theta) * pow(1 + pow(rho, 2), -0.5);
  Drho2(1, 1) -= cos(theta) * pow(rho, 2) * pow(1 + pow(rho, 2), -1.5);
  return Drho2;
}
