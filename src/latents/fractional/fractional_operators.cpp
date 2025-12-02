// fractional_operators.cpp
// C++ implementation mirroring the core logic of fractional.operators in rSPDE
//
// Inputs:
//  - C, Ci, G: FEM matrices (Eigen::SparseMatrix<double>)
//  - beta: fractional power (double)
//  - m: order of the rational approximation (int)
//  - tau: vector of tau values (Eigen::VectorXd). Can be length 1 (constant) or
//  length n.
//  - theta_kappa: parameter vector (Eigen::VectorXd)
//  - B_kappa: design matrix for log-kappa (dense Eigen::MatrixXd), so
//    kappa(s) = exp(B_kappa * theta_kappa)
//
// Outputs:
//  - Pl, Pr: assembled sparse operators as Eigen::SparseMatrix<double>
//
// Notes:
//  - This function follows the same scaling strategy as fractional.operators:
//      L <- (G + C * diag(kappa^2)) / c, with c = min(kappa)^2, and A := Ci *
//      L.
//    The final Pl is multiplied by c^beta so that the scaling cancels
//    analytically while improving numerical conditioning during assembly.
//  - For integer beta, no rational roots are needed.
//  - For fractional beta, the construction requires rational approximation
//  roots
//    rb (length m+1), rc (length m) and a scaling factor. These are provided in
//    rSPDE via internal tables (sysdata.rda) and accessed through get.roots().
//    In this standalone C++ implementation, users must supply those roots via
//    the optional overload or integrate a table-backed get_roots(). See
//    get_roots() stub below.

#include "fractional_operators.hpp"
#include "roots_tables.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

namespace rspde_cpp {

namespace {

enum class InterpType { Linear, Spline };

struct TableView {
  const double *beta;
  const double *factor;
  std::array<const double *, 4> rc_cols;
  std::array<const double *, 5> rb_cols;
  int rc_count;
  int rb_count;
  std::size_t size;
  double beta_min;
  double beta_max;
};

const TableView &table_for_order(int order) {
  using namespace detail;
  static const TableView table1{m1_beta,
                                m1_factor,
                                {m1_rc, nullptr, nullptr, nullptr},
                                {m1_rb_1, m1_rb_2, nullptr, nullptr, nullptr},
                                1,
                                2,
                                m1_size,
                                m1_beta[0],
                                m1_beta[m1_size - 1]};
  static const TableView table2{m2_beta,
                                m2_factor,
                                {m2_rc_1, m2_rc_2, nullptr, nullptr},
                                {m2_rb_1, m2_rb_2, m2_rb_3, nullptr, nullptr},
                                2,
                                3,
                                m2_size,
                                m2_beta[0],
                                m2_beta[m2_size - 1]};
  static const TableView table3{m3_beta,
                                m3_factor,
                                {m3_rc_1, m3_rc_2, m3_rc_3, nullptr},
                                {m3_rb_1, m3_rb_2, m3_rb_3, m3_rb_4, nullptr},
                                3,
                                4,
                                m3_size,
                                m3_beta[0],
                                m3_beta[m3_size - 1]};
  static const TableView table4{m4_beta,
                                m4_factor,
                                {m4_rc_1, m4_rc_2, m4_rc_3, m4_rc_4},
                                {m4_rb_1, m4_rb_2, m4_rb_3, m4_rb_4, m4_rb_5},
                                4,
                                5,
                                m4_size,
                                m4_beta[0],
                                m4_beta[m4_size - 1]};

  switch (order) {
  case 1:
    return table1;
  case 2:
    return table2;
  case 3:
    return table3;
  case 4:
    return table4;
  default:
    throw std::invalid_argument("order must be between 1 and 4");
  }
}

InterpType parse_interp(const std::string &type) {
  if (type == "linear") {
    return InterpType::Linear;
  }
  if (type == "spline") {
    return InterpType::Spline;
  }
  throw std::invalid_argument(
      "type_interp must be either 'linear' or 'spline'");
}

double clamp_beta(const TableView &tbl, double beta) {
  if (beta < tbl.beta_min) {
    return tbl.beta_min;
  }
  if (beta > tbl.beta_max) {
    return tbl.beta_max;
  }
  return beta;
}

double linear_interp(const double *x, const double *y, std::size_t n,
                     double xout) {
  if (n == 0) {
    throw std::invalid_argument("table has zero length");
  }
  if (n == 1) {
    return y[0];
  }
  if (xout <= x[0]) {
    return y[0];
  }
  if (xout >= x[n - 1]) {
    return y[n - 1];
  }
  auto it = std::lower_bound(x, x + n, xout);
  std::size_t idx = static_cast<std::size_t>(it - x);
  if (idx == 0) {
    return y[0];
  }
  if (idx >= n) {
    return y[n - 1];
  }
  if (x[idx] == xout) {
    return y[idx];
  }
  std::size_t i0 = idx - 1;
  double t = (xout - x[i0]) / (x[idx] - x[i0]);
  return y[i0] + t * (y[idx] - y[i0]);
}

struct CubicCoeffs {
  std::vector<double> b;
  std::vector<double> c;
  std::vector<double> d;
};

CubicCoeffs compute_fmm_coeffs(const double *x, const double *y,
                               std::size_t n) {
  if (n < 2) {
    throw std::invalid_argument("need at least two support points");
  }

  CubicCoeffs coeffs;
  coeffs.b.assign(n, 0.0);
  coeffs.c.assign(n, 0.0);
  coeffs.d.assign(n, 0.0);

  if (n == 2) {
    double slope = (y[1] - y[0]) / (x[1] - x[0]);
    coeffs.b[0] = coeffs.b[1] = slope;
    return coeffs;
  }

  std::vector<double> diag(n, 0.0);
  std::vector<double> offdiag(n, 0.0);
  std::vector<double> rhs(n, 0.0);

  offdiag[0] = x[1] - x[0];
  rhs[1] = (y[1] - y[0]) / offdiag[0];
  diag[0] = -offdiag[0];

  const std::size_t nm1 = n - 1;
  for (std::size_t i = 1; i < nm1; ++i) {
    offdiag[i] = x[i + 1] - x[i];
    diag[i] = 2.0 * (offdiag[i - 1] + offdiag[i]);
    rhs[i + 1] = (y[i + 1] - y[i]) / offdiag[i];
    rhs[i] = rhs[i + 1] - rhs[i];
  }

  diag[nm1] = -offdiag[nm1 - 1];
  rhs[0] = 0.0;
  rhs[nm1] = 0.0;

  if (n > 3) {
    rhs[0] = rhs[2] / (x[3] - x[1]) - rhs[1] / (x[2] - x[0]);
    rhs[nm1] = rhs[nm1 - 1] / (x[nm1] - x[nm1 - 2]) -
               rhs[nm1 - 2] / (x[nm1 - 1] - x[nm1 - 3]);
    rhs[0] = rhs[0] * offdiag[0] * offdiag[0] / (x[3] - x[0]);
    rhs[nm1] =
        -rhs[nm1] * offdiag[nm1 - 1] * offdiag[nm1 - 1] / (x[nm1] - x[nm1 - 3]);
  }

  for (std::size_t i = 1; i < n; ++i) {
    double t = offdiag[i - 1] / diag[i - 1];
    diag[i] -= t * offdiag[i - 1];
    rhs[i] -= t * rhs[i - 1];
  }

  rhs[nm1] /= diag[nm1];
  for (std::size_t i = nm1; i-- > 0;) {
    rhs[i] = (rhs[i] - offdiag[i] * rhs[i + 1]) / diag[i];
  }

  coeffs.b[nm1] = (y[nm1] - y[nm1 - 1]) / offdiag[nm1 - 1] +
                  offdiag[nm1 - 1] * (rhs[nm1 - 1] + 2.0 * rhs[nm1]);
  for (std::size_t i = 0; i < nm1; ++i) {
    coeffs.b[i] = (y[i + 1] - y[i]) / offdiag[i] -
                  offdiag[i] * (rhs[i + 1] + 2.0 * rhs[i]);
    coeffs.d[i] = (rhs[i + 1] - rhs[i]) / offdiag[i];
    coeffs.c[i] = 3.0 * rhs[i];
  }
  coeffs.c[nm1] = 3.0 * rhs[nm1];
  coeffs.d[nm1] = coeffs.d[nm1 - 1];
  return coeffs;
}

double eval_cubic(const double *x, const double *y,
                  const std::vector<double> &b, const std::vector<double> &c,
                  const std::vector<double> &d, std::size_t n, double xout) {
  if (n == 0) {
    throw std::invalid_argument("table has zero length");
  }
  if (n == 1) {
    return y[0];
  }
  if (xout <= x[0]) {
    return y[0];
  }
  if (xout >= x[n - 1]) {
    return y[n - 1];
  }
  std::size_t i = 0;
  std::size_t j = n - 1;
  while (j > i + 1) {
    std::size_t k = (i + j) / 2;
    if (xout < x[k]) {
      j = k;
    } else {
      i = k;
    }
  }
  double dx = xout - x[i];
  double tmp = d[i];
  return y[i] + dx * (b[i] + dx * (c[i] + dx * tmp));
}

double spline_interp(const double *x, const double *y, std::size_t n,
                     double xout) {
  CubicCoeffs coeffs = compute_fmm_coeffs(x, y, n);
  return eval_cubic(x, y, coeffs.b, coeffs.c, coeffs.d, n, xout);
}

double interpolate_column(const TableView &tbl, const double *column,
                          double beta, InterpType mode) {
  switch (mode) {
  case InterpType::Linear:
    return linear_interp(tbl.beta, column, tbl.size, beta);
  case InterpType::Spline:
    return spline_interp(tbl.beta, column, tbl.size, beta);
  }
  throw std::logic_error("unreachable interpolation mode");
}

} // namespace

Roots get_roots(int m, double beta, const std::string &type_interp) {
  if (m < 1 || m > 4) {
    throw std::invalid_argument("order must be one of 1, 2, 3, or 4");
  }
  double beta_adj = beta;
  if (beta_adj > 2.0) {
    beta_adj -= std::floor(beta_adj - 1.0);
  }
  const TableView &tbl = table_for_order(m);
  InterpType mode = parse_interp(type_interp);
  double beta_eval = clamp_beta(tbl, beta_adj);

  Roots r;
  r.factor = interpolate_column(tbl, tbl.factor, beta_eval, mode);
  r.rb.resize(static_cast<std::size_t>(tbl.rb_count));
  for (int i = 0; i < tbl.rb_count; ++i) {
    r.rb[static_cast<std::size_t>(i)] = interpolate_column(
        tbl, tbl.rb_cols[static_cast<std::size_t>(i)], beta_eval, mode);
  }
  r.rc.resize(static_cast<std::size_t>(tbl.rc_count));
  for (int i = 0; i < tbl.rc_count; ++i) {
    r.rc[static_cast<std::size_t>(i)] = interpolate_column(
        tbl, tbl.rc_cols[static_cast<std::size_t>(i)], beta_eval, mode);
  }
  return r;
}

// Create a sparse identity matrix of size n
inline SpMat speye(Index n) {
  SpMat I(n, n);
  I.setIdentity();
  return I;
}

// Create a sparse diagonal matrix from a dense vector
inline SpMat spdiag(const VectorXd &d) {
  // Eigen::DiagonalMatrix<double, Eigen::Dynamic> D = d.asDiagonal();
  // Convert to sparse
  SpMat S(d.size(), d.size());
  S.setZero();
  for (int i = 0; i < d.size(); ++i)
    S.insert(i, i) = d[i];
  S.makeCompressed();
  return S;
}

// Multiply a sparse matrix by a scalar (in place)
inline void smat_scale_inplace(SpMat &M, double s) { M *= s; }

// Compute product prod_i (I - alpha[i] * A) for a sparse A.
// Returns a new sparse matrix. Order is left-to-right as in R code: start with
// (I - a1 A), then multiply by (I - a2 A), etc.
inline SpMat product_linear_factors(const SpMat &A,
                                    const std::vector<double> &alpha) {
  const Index n = A.rows();
  SpMat P = speye(n) - (A * alpha.front());
  for (size_t i = 1; i < alpha.size(); ++i) {
    SpMat F = speye(n) - (A * alpha[i]);
    P = P * F; // sparse-sparse multiplication
  }
  return P;
}

// Compute A^k (k >= 0) for sparse A by repeated multiplication.
inline SpMat smat_power(const SpMat &A, int k) {
  const Index n = A.rows();
  if (k <= 0) {
    return speye(n);
  }
  SpMat P = A;
  for (int i = 1; i < k; ++i) {
    P = P * A;
  }
  return P;
}

// Build Phi = diag(1/tau). tau can be length 1 (constant) or length n.
inline SpMat build_Phi(Index n, const VectorXd &tau) {
  if (tau.size() == 1) {
    double invtau = 1.0 / tau[0];
    return spdiag(VectorXd::Constant(n, invtau));
  } else if (tau.size() == n) {
    VectorXd invtau = tau.cwiseInverse();
    return spdiag(invtau);
  } else {
    throw std::invalid_argument("tau must be length 1 or length n");
  }
}

// Core implementation, fractional beta requires rb/rc/factor.
std::pair<SpMat, SpMat> compute_fractional_operators_with_roots(
    const SpMat &C_in, const SpMat &Ci_in, const SpMat &G, double beta, int m,
    const VectorXd &tau, const VectorXd &theta_kappa, const MatrixXd &B_kappa,
    const std::vector<double> &rb, // size m+1 (only used if beta not integer)
    const std::vector<double> &rc, // size m   (only used if beta not integer)
    double roots_factor            // (only used if beta not integer)
) {
  if (m < 0) {
    throw std::invalid_argument("m must be a non-negative integer");
  }
  const Index n = C_in.rows();
  if (C_in.cols() != n || Ci_in.rows() != n || Ci_in.cols() != n ||
      G.rows() != n || G.cols() != n || B_kappa.rows() != n) {
    throw std::invalid_argument(
        "Dimension mismatch in inputs (C, Ci, G, B_kappa must be n x n (except "
        "B_kappa which is n x p))");
  }

  SpMat C_lump = C_in;
  SpMat Ci = Ci_in;
  if (C_lump.cols() != C_lump.rows()) {
    throw std::invalid_argument("C must be square");
  }
  if (Ci.cols() != Ci.rows()) {
    throw std::invalid_argument("Ci must be square");
  }
  VectorXd lump_diag = C_lump.diagonal();
  if ((lump_diag.array() <= 0).any()) {
    throw std::runtime_error("Non-positive entries in lumped C diagonal");
  }

  // kappa(s) = exp(B_kappa * theta_kappa)
  if (B_kappa.cols() != theta_kappa.size()) {
    throw std::invalid_argument("B_kappa cols must match length(theta_kappa)");
  }
  VectorXd log_kappa = B_kappa * theta_kappa;
  VectorXd kappa = log_kappa.array().exp().matrix();
  VectorXd kappa2 = kappa.array().square().matrix();

  // L = G + C_in * diag(kappa^2)  (use original C_in, matching R path)
  SpMat Kp = spdiag(kappa2);
  SpMat L = G + (C_in * Kp);

  // scale.factor = min(kappa)^2
  double scale_factor = kappa2.minCoeff();
  if (!(scale_factor > 0)) {
    throw std::runtime_error("scale_factor must be positive");
  }

  // L_scaled and A = Ci * L_scaled
  SpMat L_scaled = L;
  smat_scale_inplace(L_scaled, 1.0 / scale_factor);
  SpMat A = Ci * L_scaled;

  // Identity and Phi
  SpMat I = speye(n);
  SpMat Phi = build_Phi(n, tau);

  // Check if beta is (near) integer
  double beta_round = std::round(beta);
  bool beta_is_int = std::abs(beta - beta_round) < 1e-12;

  SpMat Pl, Pr;

  if (beta_is_int) {
    // Integer beta: Pr = I; Pl = L_scaled * A^(beta-1); then scale by
    // scale_factor^beta
    int ibeta = static_cast<int>(beta_round);
    Pr = I; // no Phi multiplication in integer case in R until the very end
            // (Phi applies regardless)
    if (ibeta <= 0) {
      // beta should be positive, but guard anyway
      Pl = I; // degenerate
    } else {
      if (ibeta == 1) {
        Pl = L_scaled;
      } else {
        SpMat Apow = smat_power(A, ibeta - 1);
        Pl = L_scaled * Apow;
      }
      // Apply final scaling: Pl <- Pl * scale_factor^beta
      double s = std::pow(scale_factor, beta);
      smat_scale_inplace(Pl, s);
    }
    // Pr <- Phi %*% Pr
    Pr = Phi * Pr;

  } else {
    // Fractional beta: need roots
    if (m == 0) {
      // No rational factors, fall back to integer-like path with m_beta only
      int m_beta = std::max(1, static_cast<int>(std::floor(beta)));
      SpMat Lp =
          C_lump *
          smat_power(A, m_beta - 1); // when m_beta==1, this is just C_lump
      Pl = Lp;                       // no (I - rb*A) factors if m==0
      // scaling
      double s = std::pow(scale_factor, beta);
      smat_scale_inplace(Pl, s);
      Pr = Phi * speye(n);
    } else {
      if (static_cast<int>(rb.size()) != m + 1 ||
          static_cast<int>(rc.size()) != m) {
        throw std::invalid_argument("rb and rc must have sizes m+1 and m, "
                                    "respectively, for fractional beta");
      }
      int m_beta = std::max(1, static_cast<int>(std::floor(beta)));

      // Pl_raw = prod_i (I - rb[i] * A)
      SpMat Pl_raw = product_linear_factors(A, rb);
      // Lp = C * A^(m_beta-1) (if m_beta==1, that power is I and Lp=C)
      SpMat Lp = C_lump * smat_power(A, m_beta - 1);
      // Pl = (scale_factor^beta / factor) * Lp * Pl_raw
      Pl = Lp * Pl_raw;
      double s = std::pow(scale_factor, beta) / roots_factor;
      smat_scale_inplace(Pl, s);

      // Pr_raw = prod_i (I - rc[i] * A)
      SpMat Pr_raw = product_linear_factors(A, rc);
      // Pr = Phi * Pr_raw
      Pr = Phi * Pr_raw;
    }
  }

  // debug
  // std::cout << "norm of Pl is " << Pl.norm() << std::endl;
  // std::cout << "norm of Pr is " << Pr.norm() << std::endl;
  return {Pl, Pr};
}

// Convenience wrapper that tries to use get_roots() when beta is fractional.
std::pair<SpMat, SpMat> compute_fractional_operators(
    const SpMat &C_in, const SpMat &Ci_in, const SpMat &G, double beta, int m,
    const VectorXd &tau, const VectorXd &theta_kappa, const MatrixXd &B_kappa) {

  // Determine if beta is fractional
  double beta_round = std::round(beta);
  bool beta_is_int = std::abs(beta - beta_round) < 1e-12;

  if (beta_is_int || m == 0) {
    // No roots needed
    static const std::vector<double> empty;
    return compute_fractional_operators_with_roots(
        C_in, Ci_in, G, beta, m, tau, theta_kappa, B_kappa, empty, empty, 1.0);
  } else {
    // Fetch roots via stub (must be implemented with embedded tables)
    Roots r = get_roots(m, beta, "linear");
    return compute_fractional_operators_with_roots(C_in, Ci_in, G, beta, m, tau,
                                                   theta_kappa, B_kappa, r.rb,
                                                   r.rc, r.factor);
  }
}

} // namespace rspde_cpp

/*
Example usage:

  using namespace rspde_cpp;
  SpMat C, Ci, G, Bk; // fill
  VectorXd tau, theta;
  double beta = 0.9;
  int m = 2;

  auto [Pl, Pr] = compute_fractional_operators(C, Ci, G, beta, m, tau, theta,
Bk);

For fractional beta, implement get_roots(m,beta) using the internal tables
from rSPDE (util.R:get.roots). Alternatively, call the ..._with_roots overload
and provide rb, rc, and factor explicitly.
*/
