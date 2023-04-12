#include "../latent.h"

Matern_2d::Matern_2d(const Rcpp::List& model_list, unsigned long seed)
  : Latent(model_list, seed),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
    alpha       (2),
    Cdiag       (C.diagonal())
  {

    // update Drho

    // update Dtheta
}

SparseMatrix<double> Matern_2d::getK(const VectorXd& alpha) const {
  Matrix2d D;
  double theta = theta_K(0);
  double rho = theta_K(1);
  double kappa1 = theta_K(2);
  double kappa2 = theta_K(3);

  D(0,0) = cos(theta) + rho*sin(theta);
  D(0,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
  D(1,0) = sin(theta) - rho*cos(theta);
  D(1,1) = cos(theta)*pow(1+pow(rho,2),0.5);

  MatrixXd Drho(2,2);
  Drho(0,0) = sin(theta);
  Drho(0,1) = -sin(theta)*rho*pow(1+pow(rho,2),-0.5);
  Drho(1,0) = -cos(theta);
  Drho(1,1) = cos(theta)*rho*pow(1+pow(rho,2),-0.5);

  MatrixXd Drho2(2,2);
  Drho2(0,0) = 0;
  Drho2(0,1) = -sin(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(0,1)+= sin(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
  Drho2(1,0) = 0;
  Drho2(1,1) = cos(theta)*pow(1+pow(rho,2),-0.5);
  Drho2(1,1) -= cos(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);

  int estimate_theta = 1;
  MatrixXd Dtheta(2,2);
  MatrixXd Dtheta2(2,2);
  if(estimate_theta == 1){
    Dtheta(0,0) = -sin(theta) + rho*cos(theta);
    Dtheta(0,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
    Dtheta(1,0) = cos(theta) + rho*sin(theta);
    Dtheta(1,1) = -sin(theta)*pow(1+pow(rho,2),0.5);

    Dtheta2(0,0) = -cos(theta) - rho*sin(theta);
    Dtheta2(0,1) = sin(theta)*pow(1+pow(rho,2),0.5);
    Dtheta2(1,0) = -sin(theta) + rho*cos(theta);
    Dtheta2(1,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
  }

  SparseMatrix<double,0,int> B,K1,K2;
  // alpha=2 case
  K1 = G + kappa1*kappa1*C;
  K2 = G + kappa2*kappa2*C;
}
SparseMatrix<double> Matern_2d::get_dK(int index, const VectorXd& alpha) const {


}

VectorXd Matern_2d::grad_theta_K() {
  VectorXd grad_theta_K = VectorXd::Zero(4);

  return grad_theta_K;
}

void Matern_2d::update_each_iter() {

}