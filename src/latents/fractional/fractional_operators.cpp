// fractional_operators.cpp
// C++ implementation mirroring the core logic of fractional.operators in rSPDE
//
// Inputs:
//  - C, Ci, G: FEM matrices (Eigen::SparseMatrix<double>)
//  - beta: fractional power (double)
//  - m: order of the rational approximation (int)
//  - tau: vector of tau values (Eigen::VectorXd). Can be length 1 (constant) or length n.
//  - theta_kappa: parameter vector (Eigen::VectorXd)
//  - B_kappa: design matrix for log-kappa (Eigen::SparseMatrix<double>), so kappa(s) = exp(B_kappa * theta_kappa)
//
// Outputs:
//  - Pl, Pr: assembled sparse operators as Eigen::SparseMatrix<double>
//
// Notes:
//  - This function follows the same scaling strategy as fractional.operators:
//      L <- (G + C * diag(kappa^2)) / c, with c = min(kappa)^2, and A := Ci * L.
//    The final Pl is multiplied by c^beta so that the scaling cancels analytically
//    while improving numerical conditioning during assembly.
//  - For integer beta, no rational roots are needed.
//  - For fractional beta, the construction requires rational approximation roots
//    rb (length m+1), rc (length m) and a scaling factor. These are provided in rSPDE
//    via internal tables (sysdata.rda) and accessed through get.roots().
//    In this standalone C++ implementation, users must supply those roots via the
//    optional overload or integrate a table-backed get_roots(). See get_roots() stub below.

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <Rcpp.h>
#include "fractional_operators.hpp"

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Index;

namespace rspde_cpp {

// Attempt to fetch roots from rSPDE::get.roots at runtime. This mirrors the R code path
// and keeps ngme2 independent of embedded tables. If rSPDE is not available,
// callers should pass rb/rc explicitly, or set m=0 to avoid fractional factors.
Roots get_roots(int m, double beta) {
  // normalize beta to (1,2] as in rSPDE (beta > 2 -> beta - floor(beta - 1))
  if (beta > 2.0) {
    beta = beta - std::floor(beta - 1.0);
  }
  try {
    Rcpp::Environment ns = Rcpp::Environment::namespace_env("rSPDE");
    if (!ns.exists("get.roots")) {
      throw std::runtime_error("rSPDE::get.roots not found");
    }
    Rcpp::Function get_roots_fun = ns["get.roots"];
    // Call with positional args (order, beta); default interpolation 'linear'
    Rcpp::List out = get_roots_fun(Rcpp::wrap(m), Rcpp::wrap(beta));
    Roots r;
    r.rb = Rcpp::as<std::vector<double>>(out["rb"]);
    r.rc = Rcpp::as<std::vector<double>>(out["rc"]);
    r.factor = Rcpp::as<double>(out["factor"]);
    // Basic sanity
    if ((int)r.rb.size() != m+1 || (int)r.rc.size() != m) {
      throw std::runtime_error("rSPDE::get.roots returned unexpected sizes");
    }
    return r;
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Unable to compute fractional roots in C++: ") + e.what());
  }
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

// Row sums of a sparse matrix (as dense vector)
inline VectorXd rowSums(const SpMat &M) {
  // rowSums = M * ones
  VectorXd ones = VectorXd::Ones(M.cols());
  VectorXd r = M * ones;
  return r;
}

// Multiply a sparse matrix by a scalar (in place)
inline void smat_scale_inplace(SpMat &M, double s) {
  M *= s;
}

// Compute product prod_i (I - alpha[i] * A) for a sparse A.
// Returns a new sparse matrix. Order is left-to-right as in R code: start with (I - a1 A), then multiply by (I - a2 A), etc.
inline SpMat product_linear_factors(const SpMat &A, const std::vector<double> &alpha) {
  const Index n = A.rows();
  SpMat P = speye(n) - (A * alpha.front());
  for (size_t i = 1; i < alpha.size(); ++i) {
    SpMat F = speye(n) - (A * alpha[i]);
    P = P * F; // sparse-sparse multiplication
  }
  return P;
}

// Fast row scaling: return M_out = D * M with D = diag(d)
inline SpMat left_scale_rows(const SpMat &M, const VectorXd &d) {
  SpMat out = M; // copy structure and values
  for (int k = 0; k < out.outerSize(); ++k) {
    for (SpMat::InnerIterator it(out, k); it; ++it) {
      it.valueRef() *= d[it.row()];
    }
  }
  return out;
}

// Check if a sparse matrix is strictly diagonal
inline bool is_diagonal_sparse(const SpMat &M) {
  for (int k = 0; k < M.outerSize(); ++k) {
    for (SpMat::InnerIterator it(M, k); it; ++it) {
      if (it.row() != it.col() && it.value() != 0.0) return false;
    }
  }
  return true;
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
    const SpMat &C_in,
    const SpMat &G,
    const VectorXd &Cdiag,
    double beta,
    int m,
    const VectorXd &tau,
    const VectorXd &theta_kappa,
    const SpMat &B_kappa,
    const std::vector<double> &rb,   // size m+1 (only used if beta not integer)
    const std::vector<double> &rc,   // size m   (only used if beta not integer)
    double roots_factor               // (only used if beta not integer)
) {
  if (m < 0) {
    throw std::invalid_argument("m must be a non-negative integer");
  }
  const Index n = C_in.rows();

  // Mass lumping as in R: build lumped C and Ci from the provided C_in,
  // but DO NOT use the lumped C to assemble L (R keeps L as provided by caller)
  // VectorXd cdiag = rowSums(C_in);
  // if ((cdiag.array() <= 0).any()) {
  //   throw std::runtime_error("Non-positive entries in rowSums(C)");
  // }
  // SpMat C_lump = spdiag(cdiag);
  // SpMat Ci = spdiag(cdiag.cwiseInverse());
  SpMat C_lump = spdiag(Cdiag);
  SpMat Ci = spdiag(Cdiag.cwiseInverse());

  // kappa(s) = exp(B_kappa * theta_kappa)
  if (B_kappa.cols() != theta_kappa.size()) {
    throw std::invalid_argument("B_kappa cols must match length(theta_kappa)");
  }
  VectorXd log_kappa = B_kappa * theta_kappa;
  // Sanitize log_kappa to avoid NaN/Inf → ensures positive finite kappa
  for (int i = 0; i < log_kappa.size(); ++i) {
    if (!std::isfinite(log_kappa[i])) log_kappa[i] = 0.0; // fallback
  }
  VectorXd kappa = log_kappa.array().exp().matrix();
  // Clamp kappa to a safe positive interval
  const double kappa_min = 1e-8;
  const double kappa_max = 1e8;
  for (int i = 0; i < kappa.size(); ++i) {
    if (!std::isfinite(kappa[i]) || kappa[i] <= 0.0) kappa[i] = kappa_min;
    else if (kappa[i] < kappa_min) kappa[i] = kappa_min;
    else if (kappa[i] > kappa_max) kappa[i] = kappa_max;
  }
  VectorXd kappa2 = kappa.array().square().matrix();

  // L = G + C_in * diag(kappa^2)
  // If C_in is diagonal (mass-lumped), avoid general sparse multiplication
  SpMat L;
  if (is_diagonal_sparse(C_in)) {
    // diag addition: diag(Cdiag * kappa^2)
    VectorXd add = Cdiag.cwiseProduct(kappa2);
    SpMat Dadd = spdiag(add);
    L = G + Dadd;
  } else {
    SpMat Kp = spdiag(kappa2);
    L = G + (C_in * Kp);
  }

  // scale.factor = min(kappa)^2
  double scale_factor = kappa2.minCoeff();
  if (!std::isfinite(scale_factor) || scale_factor <= 0) {
    // Fall back to conservative positive scaling
    scale_factor = std::max(1e-16, kappa2.array().isFinite().select(kappa2.array(), 1.0).minCoeff());
  }

  // L_scaled and A = Ci * L_scaled
  SpMat L_scaled = L;
  smat_scale_inplace(L_scaled, 1.0 / scale_factor);
  // A = Ci * L_scaled; when Ci is diagonal, do fast row scaling
  SpMat A;
  if (is_diagonal_sparse(Ci)) {
    VectorXd invC = Cdiag.cwiseInverse();
    A = left_scale_rows(L_scaled, invC);
  } else {
    A = Ci * L_scaled;
  }

  // Identity and Phi
  SpMat I = speye(n);
  SpMat Phi = build_Phi(n, tau);

  // Check if beta is (near) integer
  double beta_round = std::round(beta);
  bool beta_is_int = std::abs(beta - beta_round) < 1e-12;

  SpMat Pl, Pr;

  if (beta_is_int) {
    // Integer beta: Pr = I; Pl = L_scaled * A^(beta-1); then scale by scale_factor^beta
    int ibeta = static_cast<int>(beta_round);
    Pr = I; // no Phi multiplication in integer case in R until the very end (Phi applies regardless)
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
      SpMat Lp = C_lump * smat_power(A, m_beta - 1); // when m_beta==1, this is just C_lump
      Pl = Lp; // no (I - rb*A) factors if m==0
      // scaling
      double s = std::pow(scale_factor, beta);
      smat_scale_inplace(Pl, s);
      Pr = Phi * speye(n);
    } else {
      if (static_cast<int>(rb.size()) != m + 1 || static_cast<int>(rc.size()) != m) {
        throw std::invalid_argument("rb and rc must have sizes m+1 and m, respectively, for fractional beta");
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

  return {Pl, Pr};
}

// Convenience wrapper that tries to use get_roots() when beta is fractional.
std::pair<SpMat, SpMat> compute_fractional_operators(
    const SpMat &C_in,
    const SpMat &G,
    const VectorXd &Cdiag,
    double beta,
    int m,
    const VectorXd &tau,
    const VectorXd &theta_kappa,
    const SpMat &B_kappa) {

  // Determine if beta is fractional
  double beta_round = std::round(beta);
  bool beta_is_int = std::abs(beta - beta_round) < 1e-12;

  if (beta_is_int || m == 0) {
    // No roots needed
    static const std::vector<double> empty;
    return compute_fractional_operators_with_roots(C_in, G, Cdiag, beta, m, tau, theta_kappa, B_kappa, empty, empty, 1.0);
  } else {
    // Fetch roots via stub (must be implemented with embedded tables)
    Roots r = get_roots(m, beta);
    return compute_fractional_operators_with_roots(C_in, G, Cdiag, beta, m, tau, theta_kappa, B_kappa, r.rb, r.rc, r.factor);
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

  auto [Pl, Pr] = compute_fractional_operators(C, Ci, G, beta, m, tau, theta, Bk);

For fractional beta, implement get_roots(m,beta) using the internal tables
from rSPDE (util.R:get.roots). Alternatively, call the ..._with_roots overload
and provide rb, rc, and factor explicitly.
*/
