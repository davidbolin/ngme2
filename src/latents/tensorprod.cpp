// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x

#include "../operator.h"
#include <algorithm>
#include "MatrixAlgebra.h"
#include <unsupported/Eigen/KroneckerProduct>

// tensor product for the C G class
Tensor_prod::Tensor_prod(const Rcpp::List &operator_list)
    : Operator(operator_list),
      first(OperatorFactory::create(operator_list["first"])),
      second(OperatorFactory::create(operator_list["second"])),
      n_theta_1(first->get_n_theta_K()), n_theta_2(second->get_n_theta_K()) {}

void Tensor_prod::build_KZ(const VectorXd &theta_K) {
  // report the time for this function
  double time = 0;
  auto timer_computeg = std::chrono::steady_clock::now();
  // std::cout << "update K now" << std::endl;
  first->build_KZ(theta_K.segment(0, n_theta_1));
  // std::cout << "update K now1" << std::endl;
  second->build_KZ(theta_K.segment(n_theta_1, n_theta_2));
  // std::cout << "new K size = " << first->getK().rows() *
  // second->getK().rows() << " " << first->getK().cols() *
  // second->getK().cols() << std::endl;

  // use Eigen kronecker product
  KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double>>
      kroneckerEigen(first->getK(), second->getK());

  kroneckerEigen.evalTo(K);
  // std::cout << "update K now3" << std::endl;

  time = std::chrono::duration_cast<std::chrono::milliseconds>(
             std::chrono::steady_clock::now() - timer_computeg)
             .count();
  // std::cout << "size and time for kronecker product is " << K.rows() << " "
  // << K.cols() << " " << time << std::endl; Z remains identity (set once in
  // Operator base)
}


// Traces of a tensor product operator, without ever factorizing K.
//
// With K = K_1 (x) K_2 the derivatives keep the Kronecker structure,
//     dK/dtheta_1j = (dK_1/dtheta_1j) (x) K_2,
//     dK/dtheta_2j = K_1 (x) (dK_2/dtheta_2j),
// and (A (x) B)^-1 = A^-1 (x) B^-1, (A (x) B)(C (x) D) = AC (x) BD give
//     tr(K^-1 (dK_1 (x) K_2)) = tr(K_1^-1 dK_1) tr(I_n2) = n_2 tr(K_1^-1 dK_1)
//     tr(K^-1 (K_1 (x) dK_2)) = n_1 tr(K_2^-1 dK_2).
// The same algebra closes the H_K block. For j, k in the same factor,
//     tr(K^-1 dK_k K^-1 dK_j) = n_other tr(F^-1 dF_k F^-1 dF_j)
//     tr(K^-1 d2K_jk)         = n_other tr(F^-1 d2F_jk),
// so that block is just n_other times the factor's own H_K block. For j in one
// factor and k in the other the two terms are equal and opposite,
//     tr(K^-1 dK_k K^-1 dK_j) = tr(K_1^-1 dK_1j) tr(K_2^-1 dK_2k)
//     d2K/dtheta_1j dtheta_2k = dK_1j (x) dK_2k, same trace,
// so the cross block vanishes identically.
bool Tensor_prod::compute_traces_structured(const VectorXd &theta,
                                            const UpdateOptions &opts) {
  const int n1 = static_cast<int>(first->getK().rows());
  const int n2 = static_cast<int>(second->getK().rows());
  if (n1 <= 0 || n2 <= 0 || n_theta_1 + n_theta_2 != n_theta_K)
    return false;

  auto budget = [&](int n_self, int n_other) {
    const long long want =
        static_cast<long long>(std::max(1, opts.n_trace_iter)) * n_other;
    return static_cast<int>(std::min<long long>(want, n_self));
  };

  UpdateOptions sub = opts;
  sub.compute_K = true;
  sub.compute_Z = true;
  sub.compute_dK = true;
  sub.compute_dZ = false;
  sub.compute_d2K = opts.compute_HK_trace;
  sub.compute_d2Z = false;
  sub.compute_trace = true;

  const bool have_mask =
      static_cast<int>(opts.fix_mask_thetaK.size()) == n_theta_K;

  UpdateOptions o1 = sub;
  o1.n_trace_iter = budget(n1, n2);
  if (have_mask)
    o1.fix_mask_thetaK.assign(opts.fix_mask_thetaK.begin(),
                              opts.fix_mask_thetaK.begin() + n_theta_1);
  else
    o1.fix_mask_thetaK.clear();

  UpdateOptions o2 = sub;
  o2.n_trace_iter = budget(n2, n1);
  if (have_mask)
    o2.fix_mask_thetaK.assign(opts.fix_mask_thetaK.begin() + n_theta_1,
                              opts.fix_mask_thetaK.end());
  else
    o2.fix_mask_thetaK.clear();

  // Also restores each factor's K to base theta: the numeric differencing in
  // Operator::update_all rebuilds the composite from perturbed factors and
  // restores only the composite, leaving the factors at the last perturbed
  // value.
  first->update_all(theta.segment(0, n_theta_1), o1);
  second->update_all(theta.segment(n_theta_1, n_theta_2), o2);

  if (!first->traces_ready() || !second->traces_ready())
    return false;
  const VectorXd &t1 = first->get_trace_trK();
  const VectorXd &t2 = second->get_trace_trK();
  if (t1.size() != n_theta_1 || t2.size() != n_theta_2)
    return false;

  if (trace_vals.size() != n_theta_K)
    trace_vals = VectorXd::Zero(n_theta_K);
  for (int j = 0; j < n_theta_1; ++j)
    trace_vals(j) = n2 * t1(j);
  for (int j = 0; j < n_theta_2; ++j)
    trace_vals(n_theta_1 + j) = n1 * t2(j);

  if (opts.compute_HK_trace) {
    const MatrixXd &H1 = first->get_HK_trace();
    const MatrixXd &H2 = second->get_HK_trace();
    if (H1.rows() != n_theta_1 || H2.rows() != n_theta_2)
      return false;
    if (HK_trace.rows() != n_theta_K || HK_trace.cols() != n_theta_K)
      HK_trace = MatrixXd::Zero(n_theta_K, n_theta_K);
    else
      HK_trace.setZero(); // cross blocks stay zero
    HK_trace.topLeftCorner(n_theta_1, n_theta_1) = n2 * H1;
    HK_trace.bottomRightCorner(n_theta_2, n_theta_2) = n1 * H2;
  }
  return true;
}

// Non-separable Space-time model

Spacetime::Spacetime(const Rcpp::List &operator_list)
    : Operator(operator_list),
      Ct_diag(Rcpp::as<VectorXd>(operator_list["Ct_diag"])),
      Cs_diag(Rcpp::as<VectorXd>(operator_list["Cs_diag"])),
      BtCs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["BtCs"])),
      Gs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Gs"])),
      Ct(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Ct"])),
      Cs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Cs"])),
      Bx(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Bx"])),
      By(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["By"])),
      S(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["S"])),
      Bs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Bs"])),
      Hxx(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hxx"])),
      Hyy(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hyy"])),
      Hxy(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hxy"])),
      Hyx(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hyx"])),
      // B_gamma_x (Rcpp::as<MatrixXd> (operator_list["B_gamma_x"])),
      // B_gamma_y (Rcpp::as<MatrixXd> (operator_list["B_gamma_y"])),
      B_gamma_x_list_input(
          Rcpp::as<Rcpp::List>(operator_list["B_gamma_x_list"])),
      B_gamma_y_list_input(
          Rcpp::as<Rcpp::List>(operator_list["B_gamma_y_list"])),
      theta_gamma_x(Rcpp::as<VectorXd>(operator_list["theta_gamma_x"])),
      theta_gamma_y(Rcpp::as<VectorXd>(operator_list["theta_gamma_y"])),
      n_theta_gamma_x(Rcpp::as<int>(operator_list["n_theta_gamma_x"])),
      n_theta_gamma_y(Rcpp::as<int>(operator_list["n_theta_gamma_y"])),
      lambda(Rcpp::as<double>(operator_list["lambda"])),
      alpha(Rcpp::as<double>(operator_list["alpha"])),
      method(Rcpp::as<string>(operator_list["method"])),
      stabilization(Rcpp::as<bool>(operator_list["stabilization"])),
      fix_gamma(Rcpp::as<bool>(operator_list["fix_gamma"])),
      shared_theta_gamma(Rcpp::as<bool>(operator_list["shared_theta_gamma"])),
      nt(Rcpp::as<int>(operator_list["nt"])), B_gamma_x_list(nt - 1),
      B_gamma_y_list(nt - 1) {
  // turn B_gamma_x_list_input into B_gamma_x_list
  if (!fix_gamma) {
    for (int i = 0; i < nt - 1; i++) {
      B_gamma_x_list[i] = Rcpp::as<MatrixXd>(B_gamma_x_list_input[i]);
      B_gamma_y_list[i] = Rcpp::as<MatrixXd>(B_gamma_y_list_input[i]);
    }
  }
}

void Spacetime::build_KZ(const VectorXd &theta_K) {
  K.setZero();
  double c = exp(theta_K[0]);
  double kappa = exp(theta_K[1]);

  if (!fix_gamma) {
    if (shared_theta_gamma) {
      theta_gamma_x = theta_K.segment(2, n_theta_gamma_x);
      theta_gamma_y = theta_gamma_x; // Use same parameter for both x and y
    } else {
      theta_gamma_x = theta_K.segment(2, n_theta_gamma_x);
      theta_gamma_y = theta_K.segment(2 + n_theta_gamma_x, n_theta_gamma_y);
    }

    std::vector<VectorXd> gamma_x_list(nt - 1);
    std::vector<VectorXd> gamma_y_list(nt - 1);
    for (int i = 0; i < nt - 1; i++) {
      gamma_x_list[i] = B_gamma_x_list[i] * theta_gamma_x;
      gamma_y_list[i] = B_gamma_y_list[i] * theta_gamma_y;
    }

    // Build Bs_list, S_list, Ls_list
    std::vector<SparseMatrix<double, 0, int>> Bs_list(nt - 1);
    std::vector<SparseMatrix<double, 0, int>> S_list(nt - 1);
    std::vector<SparseMatrix<double, 0, int>> Ls_list(nt - 1);
    for (int i = 0; i < nt - 1; i++) {
      Bs_list[i] =
          gamma_x_list[i].asDiagonal() * Bx * gamma_x_list[i].asDiagonal() +
          gamma_y_list[i].asDiagonal() * By * gamma_y_list[i].asDiagonal();

      Ls_list[i] = kappa * kappa * Cs + lambda * Gs + Bs_list[i];

      if (alpha == 4) {
        Ls_list[i] = Ls_list[i] * Cs_diag.cwiseInverse().asDiagonal() *
                     Ls_list[i].transpose();
      }

      if (stabilization) {
        // Compute S_list[i]
        if (gamma_x_list[i].norm() < 1e-8 && gamma_y_list[i].norm() < 1e-8) {
          S_list[i] = SparseMatrix<double, 0, int>(Ls_list[i].rows(),
                                                   Ls_list[i].cols());
          S_list[i].setZero();
        } else {
          VectorXd gamma_xx = gamma_x_list[i].array().square();
          VectorXd gamma_yy = gamma_y_list[i].array().square();
          VectorXd gamma_xy = gamma_x_list[i].array() * gamma_y_list[i].array();

          S_list[i] =
              gamma_xx.asDiagonal() * Hxx * gamma_xx.asDiagonal() +
              gamma_yy.asDiagonal() * Hyy * gamma_yy.asDiagonal() +
              gamma_xy.asDiagonal() * (Hxy + Hyx) * gamma_xy.asDiagonal();

          double gamma_norm = sqrt((gamma_x_list[i].array().square() +
                                    gamma_y_list[i].array().square())
                                       .sum());
          S_list[i] = Cs_diag.asDiagonal() * S_list[i] / gamma_norm;
        }
        // Add S_list[i] to Ls_list[i]
        Ls_list[i] = Ls_list[i] + S_list[i];
      }
    }

    for (int i = 1; i < Ct.rows(); ++i) {
      setSparseBlock(&K, i * Ls_list[i - 1].rows(), i * Ls_list[i - 1].cols(),
                     Ls_list[i - 1]);
    }
  } else { // fix_gamma = TRUE
    SparseMatrix<double> Ls = (kappa * kappa * Cs + lambda * Gs + Bs);
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
  //   KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> >
  //   kroneckerEigen(Ct, Ls);

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

  //   // KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> >
  //   kroneckerEigen(I_sparse, Ls);
  //   // K = kroneckerEigen.eval();

  //   // 2nd approach: build Bdiag(0, Ls, ..., Ls) directly
  //   K.setZero();
  //   for (int i = 0; i < nt-1; ++i) {
  //     setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls_list[i]);
  //   }
  // }
}

void Spacetime::update_dK(const VectorXd &theta_K) {}
