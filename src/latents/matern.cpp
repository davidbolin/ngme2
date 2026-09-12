/*
    Matern model (stationary or non-stationary kappa):
        kappa(s) = exp(B_theta_K * theta_K)
        alpha is the smoothness parameter
        K = kappa^2 * C + G for integer alpha (with optional fractional case)
*/

#include "../operator.h"
#include "fractional/fractional_operators.hpp"
#include <stdexcept>

Matern::Matern(const Rcpp::List &operator_list)
    : Operator(operator_list),
      G(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["G"])),
      C(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["C"])),
      Ci(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Ci"])),
      Bk_dense(Rcpp::as<MatrixXd>(operator_list["B_K"])),
      alpha(Rcpp::as<double>(operator_list["alpha"])), Cdiag(C.diagonal()) {
  // Stationarity flag and basis matrix for log kappa (dense)
  stationary = operator_list.containsElementNamed("stationary")
                   ? Rcpp::as<bool>(operator_list["stationary"])
                   : true;

  // Optional fractional controls
  fix_alpha = operator_list.containsElementNamed("fix_alpha")
                  ? Rcpp::as<bool>(operator_list["fix_alpha"])
                  : true;
  m = operator_list.containsElementNamed("rational_order")
          ? Rcpp::as<int>(operator_list["rational_order"])
          : 0;
  dim = operator_list.containsElementNamed("spatial_dim")
            ? Rcpp::as<int>(operator_list["spatial_dim"])
            : 2;

  const int n = G.rows();
  if (Bk_dense.rows() != n) {
    Rcpp::stop("B_K has %d rows but expected %d",
               static_cast<int>(Bk_dense.rows()), n);
  }
}

void Matern::build_KZ(const VectorXd &theta_K) {
  using namespace rspde_cpp;
  // theta_K layout:
  //   if (!fix_alpha): [eta_alpha, theta_K ...]
  //   if ( fix_alpha): [theta_K ...]
  int offset = 0;
  if (!fix_alpha) {
    double eta_alpha = theta_K(0);
    double L = 0.5 * static_cast<double>(dim);
    double sig = 1.0 / (1.0 + std::exp(-eta_alpha));
    alpha = L + (4.0 - L) * sig; // (L,4)
    offset = 1;
  }

  const int n_kappa = static_cast<int>(theta_K.size()) - offset;
  if (n_kappa <= 0) {
    throw std::invalid_argument("theta_K must contain kappa coefficients");
  }
  if (Bk_dense.cols() != n_kappa) {
    Rcpp::stop("Length of theta_K (%d) does not match B_K columns (%d)",
               n_kappa, static_cast<int>(Bk_dense.cols()));
  }
  VectorXd theta_kappa = theta_K.segment(offset, n_kappa);

  // kappa(s) = exp(B_theta_K * theta_kappa)
  VectorXd log_kappa = Bk_dense * theta_kappa;
  VectorXd kappa = log_kappa.array().exp();
  VectorXd kappa2 = kappa.array().square();

  if (std::abs(alpha - 2) < 1e-6) {
    // Integer case alpha = 2
    if (stationary) {
      SparseMatrix<double> KCK = (C * kappa2.asDiagonal()).eval();
      K = (G + KCK);
    } else {
      SparseMatrix<double> Dk(kappa.size(), kappa.size());
      Dk = kappa.asDiagonal();
      K = (Dk * C * Dk + G).eval();
    }
    Z.setIdentity();
  } else if (std::abs(alpha - 4) < 1e-6) {
    // Integer case alpha = 4
    if (stationary) {
      SparseMatrix<double> KCK = (C * kappa2.asDiagonal()).eval();
      K = (G + KCK) * Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
    } else {
      SparseMatrix<double> Dk(kappa.size(), kappa.size());
      Dk = kappa.asDiagonal();
      SparseMatrix<double> L = (Dk * C * Dk + G).eval();
      K = (L * Cdiag.cwiseInverse().asDiagonal() * L).eval();
    }
    Z.setIdentity();
  } else {
    // Fractional case
    double beta = alpha / 2.0;
    // tau defaults to 1 for now
    VectorXd tau(1);
    tau(0) = 1.0;
    auto pairKZ = compute_fractional_operators(C, Ci, G, beta, m, tau,
                                               theta_kappa, Bk_dense);
    K = pairKZ.first;  // Pl
    Z = pairKZ.second; // Pr
  }
}

// ---------------------------------------------------------------------------
// Closed-form derivatives of the Matern operator, integer alpha.
//
//   kappa = exp(B_kappa theta),   A = G + C diag(kappa^2)   [stationary]
//                                 A = G + Dk C Dk           [non-stationary]
//   alpha = 2:  K = A
//   alpha = 4:  K = A Cdiag^-1 A
//
// theta enters only through the diagonal kappa, and
//   d kappa^2_i / d theta_j    = 2 kappa^2_i B_ij
//   d2 kappa^2_i / dtheta_j dtheta_k = 4 kappa^2_i B_ij B_ik,
// so a derivative of A is one diagonal rescaling and a derivative of the
// alpha = 4 operator is two sparse products.
// ---------------------------------------------------------------------------

namespace {
// Eigen assigns a DiagonalWrapper to a SparseMatrix but will not construct one
// from it, so the two-step form is wrapped once here rather than repeated.
inline SparseMatrix<double> spdiag(const VectorXd &v) {
  SparseMatrix<double> D(v.size(), v.size());
  D = v.asDiagonal();
  return D;
}
} // namespace

bool Matern::matern_dA(const VectorXd &theta_K, SparseMatrix<double> &A,
                       std::vector<SparseMatrix<double>> &dA,
                       std::vector<std::vector<SparseMatrix<double>>> *d2A)
    const {
  if (!fix_alpha)
    return false; // theta_K(0) moves alpha; not handled here
  const bool a2 = std::abs(alpha - 2) < 1e-6;
  const bool a4 = std::abs(alpha - 4) < 1e-6;
  if (!a2 && !a4)
    return false; // fractional

  const int p = (int)theta_K.size();
  if (Bk_dense.cols() != p || p <= 0)
    return false;

  const VectorXd kappa = (Bk_dense * theta_K).array().exp();
  const VectorXd kappa2 = kappa.array().square();

  // P is the kappa-dependent part of A; G is constant, so dA = dP.
  //
  // Every derivative of P is an elementwise rescaling of P, with no matrix
  // products at all. Entry (i,l) is
  //   stationary      P_il = C_il kappa_l^2        -> d/dtheta_j = P_il 2 B_lj
  //   non-stationary  P_il = kappa_i C_il kappa_l  -> d/dtheta_j = P_il (B_ij + B_lj)
  // and in both cases the second derivative is P_il w_j w_k with the same
  // per-entry weight w. So one pass over nnz(P) per parameter gives every
  // first derivative and one multiply per entry gives every second.
  SparseMatrix<double> P;
  if (stationary) {
    P = C * kappa2.asDiagonal();
  } else {
    const SparseMatrix<double> Dk = spdiag(kappa);
    P = Dk * C * Dk;
  }
  P.makeCompressed();
  A = G + P;

  const int nz = (int)P.nonZeros();
  const int *inner = P.innerIndexPtr();
  const int *outer = P.outerIndexPtr();

  // w[j][t] = weight of stored entry t for parameter j.
  std::vector<std::vector<double>> w(p, std::vector<double>(nz, 0.0));
  for (int c = 0; c < P.outerSize(); ++c)
    for (int t = outer[c]; t < outer[c + 1]; ++t) {
      const int r = inner[t];
      for (int j = 0; j < p; ++j)
        w[j][t] = stationary ? 2.0 * Bk_dense(c, j)
                             : Bk_dense(r, j) + Bk_dense(c, j);
    }

  dA.assign(p, SparseMatrix<double>());
  for (int j = 0; j < p; ++j) {
    dA[j] = P; // same pattern, values rescaled below
    double *v = dA[j].valuePtr();
    const double *pv = P.valuePtr();
    for (int t = 0; t < nz; ++t)
      v[t] = pv[t] * w[j][t];
  }

  if (d2A) {
    d2A->assign(p, std::vector<SparseMatrix<double>>(p));
    for (int j = 0; j < p; ++j)
      for (int k = j; k < p; ++k) {
        SparseMatrix<double> out = P;
        double *v = out.valuePtr();
        const double *pv = P.valuePtr();
        for (int t = 0; t < nz; ++t)
          v[t] = pv[t] * w[j][t] * w[k][t];
        (*d2A)[j][k] = out;
        if (k != j)
          (*d2A)[k][j] = out;
      }
  }
  return true;
}

bool Matern::update_dKdZ(const VectorXd &theta_K) {
  static const bool disabled = [] {
    const char *e = std::getenv("NGME_MATERN_NUMERIC_DK");
    return e && *e && std::string(e) != "0";
  }();
  if (disabled)
    return false;

  SparseMatrix<double> A;
  std::vector<SparseMatrix<double>> dA;
  if (!matern_dA(theta_K, A, dA, nullptr))
    return false;
  const int p = (int)theta_K.size();

  if ((int)dK.size() != p)
    dK.assign(p, SparseMatrix<double>(K.rows(), K.cols()));
  if ((int)dZ.size() != p)
    dZ.assign(p, SparseMatrix<double>(K.rows(), K.cols()));

  const bool a4 = std::abs(alpha - 4) < 1e-6;
  for (int j = 0; j < p; ++j) {
    if (a4) {
      // K = A D^-1 A  ->  dK = dA D^-1 A + A D^-1 dA
      const auto Di = Cdiag.cwiseInverse().asDiagonal();
      dK[j] = SparseMatrix<double>(dA[j] * Di * A) +
              SparseMatrix<double>(A * Di * dA[j]);
    } else {
      dK[j] = dA[j];
    }
    dK[j].makeCompressed();
    // Z is the identity for both integer cases, so it carries no parameters.
    dZ[j].setZero();
  }
  return true;
}

bool Matern::update_d2Kd2Z(const VectorXd &theta_K) {
  static const bool disabled = [] {
    const char *e = std::getenv("NGME_MATERN_NUMERIC_DK");
    return e && *e && std::string(e) != "0";
  }();
  if (disabled)
    return false;

  SparseMatrix<double> A;
  std::vector<SparseMatrix<double>> dA;
  std::vector<std::vector<SparseMatrix<double>>> d2A;
  if (!matern_dA(theta_K, A, dA, &d2A))
    return false;
  const int p = (int)theta_K.size();

  if ((int)d2K.size() != p)
    d2K.assign(p, std::vector<SparseMatrix<double>>(
                      p, SparseMatrix<double>(K.rows(), K.cols())));
  if ((int)d2Z.size() != p)
    d2Z.assign(p, std::vector<SparseMatrix<double, 0, int>>(
                      p, SparseMatrix<double, 0, int>(K.rows(), K.cols())));

  const bool a4 = std::abs(alpha - 4) < 1e-6;
  for (int j = 0; j < p; ++j)
    for (int k = j; k < p; ++k) {
      SparseMatrix<double> out;
      if (a4) {
        // d2(A D^-1 A) = d2A D^-1 A + dA_j D^-1 dA_k + dA_k D^-1 dA_j
        //                + A D^-1 d2A
        const auto Di = Cdiag.cwiseInverse().asDiagonal();
        out = SparseMatrix<double>(d2A[j][k] * Di * A) +
              SparseMatrix<double>(dA[j] * Di * dA[k]) +
              SparseMatrix<double>(dA[k] * Di * dA[j]) +
              SparseMatrix<double>(A * Di * d2A[j][k]);
      } else {
        out = d2A[j][k];
      }
      out.makeCompressed();
      d2K[j][k] = out;
      d2Z[j][k].setZero();
      if (k != j) {
        d2K[k][j] = out;
        d2Z[k][j].setZero();
      }
    }
  return true;
}
