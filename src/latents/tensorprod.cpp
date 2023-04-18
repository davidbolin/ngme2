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
  // std::cout << "start init" << std::endl;

  unsigned long latent_seed = latent_rng();
  Rcpp::List group = model_list["left"];
  left = LatentFactory::create(group, latent_seed);
  Rcpp::List model_right = model_list["right"];
  right = LatentFactory::create(model_right, latent_seed);

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

  // std::cout << "Finish init" << std::endl;
}

VectorXd Tensor_prod::grad_theta_K() {
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

  for (int i = 0; i < n_rep; i++) {
    VectorXd W = Ws[i];
    VectorXd prevW = prevWs[i];
    VectorXd V = vars[i].getV();
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

    for (int j = 0; j < n1 + n2; j++) {
      SparseMatrix<double> dK;
      if (j < n1) {
        SparseMatrix<double> dKl = left->get_dK_by_index(j);
        dK = kronecker(dKl, K2);
      } else {
        SparseMatrix<double> dKr = right->get_dK_by_index(j - n1);
        dK = kronecker(K1, dKr);
      }
      // compute trace
      double trace = 0;
      if (!symmetricK) {
        lu_solver_K.computeKTK(K);
        trace = lu_solver_K.trace(dK);
      } else {
        chol_solver_K.compute(K);
        trace = chol_solver_K.trace(dK);
      }
      // compute rhs
      double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
      grad(j) += trace - tmp;
    }
  }

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
  VectorXd theta_K_l = theta_K.segment(0, n1);
  VectorXd theta_K_r = theta_K.segment(n1, n2);

  if (index < n1) {
    SparseMatrix<double> dk_l = left->get_dK(index, theta_K_l);
    SparseMatrix<double> K_r = right->getK(theta_K_r);
    return kronecker(dk_l, K_r);
  } else {
    SparseMatrix<double> dk_r = right->get_dK(index - n1, theta_K_r);
    SparseMatrix<double> K_l = left->getK(theta_K_l);
    return kronecker(K_l, dk_r);
  }
}

// // after set parameters, update K and K_rep
// void Tensor_prod::update_each_iter() {
//   int n1 = left->get_n_theta_K();
//   int n2 = right->get_n_theta_K();
//   // setK for left and right
// // std::cout << "Theta K = " << theta_K << std::endl;
//   // wrong here
//   // left->set_theta_K(theta_K.segment(0, n1));
//   // right->set_theta_K(theta_K.segment(n1, n2));

//   SparseMatrix<double> K1 = left->getK(theta_K.segment(0, n1));
//   SparseMatrix<double> K2 = right->getK(theta_K.segment(n1, n2));
//   K = kronecker(K1, K2);
//   for (int i=0; i < n_rep; i++)
//     setSparseBlock(&K_rep, i*V_size, i*W_size, K);

//   // SparseMatrix<double> Q = K.transpose() * K;

//   // dKtx_t = dKt_t %x% Kx
//   // dKtx_x = kt %x% dKx_x

//   // for now, assume 1 parameter for each only
//   // dK  = kronecker(left->get_dK(0), K2);
//   // dK2 = kronecker(K1, right->get_dK(0));
// }

// further, latent.get_trace(index) should be implemented


// ------ iid model -----

Iid::Iid(const Rcpp::List& model_list, unsigned long seed)
  : Latent(model_list, seed),
    I(Rcpp::as< SparseMatrix<double,0,int> > (model_list["K"])) {
    K = I;
}

SparseMatrix<double> Iid::getK(const VectorXd& alpha) const {return I;};
SparseMatrix<double> Iid::get_dK(int index, const VectorXd& alpha) const {return 0*I;};
VectorXd Iid::grad_theta_K() {return VectorXd::Zero(1);};
void Iid::update_each_iter() {};


// ------- bivariate model -------
Bivar::Bivar(const Rcpp::List& model_list, unsigned long seed)
  : Latent(model_list, seed) {

  Rcpp::List l1 = model_list["m1"];
  m1 = LatentFactory::create(l1, latent_rng());
  Rcpp::List l2 = model_list["m2"];
  m2 = LatentFactory::create(l2, latent_rng());
  n = m1->get_W_size();

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

  std::cout << "Finish init bivariate" << std::endl;
}

SparseMatrix<double> Bivar::getK(const VectorXd& theta_K) const {
  double theta = theta_K(0);
  double rho = theta_K(1);
  VectorXd theta_K1 = theta_K.segment(2, m1->get_n_theta_K());
  VectorXd theta_K2 = theta_K.segment(2+m1->get_n_theta_K(), m2->get_n_theta_K());

  Matrix2d D = getD(theta, rho);
  SparseMatrix<double> K1 = m1->getK(theta_K1);
  SparseMatrix<double> K2 = m2->getK(theta_K2);

  SparseMatrix<double> K (2*n, 2*n);
  SparseMatrix<double> K00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * K1;
  SparseMatrix<double> K01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * K2;
  SparseMatrix<double> K10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * K1;
  SparseMatrix<double> K11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * K2;

  setSparseBlock(&K, 0, 0, K00);
  setSparseBlock(&K, 0, n, K01);
  setSparseBlock(&K, n, 0, K10);
  setSparseBlock(&K, n, n, K11);

  return K;
}

SparseMatrix<double> Bivar::get_dK(int index, const VectorXd& theta_K) const {
  double theta = theta_K(0);
  double rho = theta_K(1);
  VectorXd theta_K1 = theta_K.segment(2, m1->get_n_theta_K());
  VectorXd theta_K2 = theta_K.segment(2+m1->get_n_theta_K(), m2->get_n_theta_K());

  SparseMatrix<double> dK (2*n, 2*n);
  if (index == 1 || index == 2) {
    // d_theta K = Dtheta * K
    SparseMatrix<double> K1 = m1->getK(theta_K1);
    SparseMatrix<double> K2 = m2->getK(theta_K2);
    Matrix2d dD;
    if (index == 1)
      dD = get_dD_theta(theta, rho);
    else
      dD = get_dD_rho(theta, rho);

    SparseMatrix<double> dK00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * K1;
    SparseMatrix<double> dK01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * K2;
    SparseMatrix<double> dK10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * K1;
    SparseMatrix<double> dK11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * K2;

    setSparseBlock(&dK, 0, 0, dK00);
    setSparseBlock(&dK, 0, n, dK01);
    setSparseBlock(&dK, n, 0, dK10);
    setSparseBlock(&dK, n, n, dK11);
  } else if (index < 2 + m1->get_n_theta_K()) {
    dK.setZero();
    SparseMatrix<double> dK1 = m1->get_dK(index-2, theta_K1);
    SparseMatrix<double> dK00 = VectorXd::Constant(n, D(0,0)).asDiagonal() * dK1;
    SparseMatrix<double> dK10 = VectorXd::Constant(n, D(1,0)).asDiagonal() * dK1;
    setSparseBlock(&dK, 0, 0, dK00);
    setSparseBlock(&dK, n, 0, dK10);
  } else {
    dK.setZero();
    SparseMatrix<double> dK2 = m2->get_dK(index-2-m1->get_n_theta_K(), theta_K2);
    SparseMatrix<double> dK01 = VectorXd::Constant(n, D(0,1)).asDiagonal() * dK2;
    SparseMatrix<double> dK11 = VectorXd::Constant(n, D(1,1)).asDiagonal() * dK2;
    setSparseBlock(&dK, 0, n, dK01);
    setSparseBlock(&dK, n, n, dK11);
  }

  return dK;
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