#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

// #include<Eigen/IterativeLinearSolvers>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>
#include <random>

#include "include/solver.h"
#include "include/timer.h"

using Eigen::Matrix2d;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using std::exp;
using std::log;
using std::pow;
using std::string;
using std::vector;

// New enums and options for unified K/Z updates and traces
enum class DiffMode { Forward = 0, Central = 1 };

struct UpdateOptions {
  bool compute_K{true};
  bool compute_Z{true};
  bool compute_dK{false};
  bool compute_dZ{false};
  bool compute_d2K{false};
  bool compute_d2Z{false};
  bool compute_HK_trace{false};
  bool prefer_analytic_dK{true};
  bool prefer_analytic_dZ{true};
  bool prefer_analytic_d2K{false};
  bool prefer_analytic_d2Z{false};
  int n_trace_iter{8};
  int solver_type{0};
  double eps_dK{1e-4};
  double eps_dZ{1e-4};
  DiffMode diff_dK_mode{DiffMode::Forward};
  DiffMode diff_dZ_mode{DiffMode::Forward};
  std::vector<bool> fix_mask_thetaK; // optional; size = n_theta_K
};

class Operator {
protected:
  VectorXd h;
  int n_theta_K;
  bool zero_trace, symmetric;
  string generic_type;

  SparseMatrix<double> K;
  vector<SparseMatrix<double>> dK;
  std::vector<std::vector<SparseMatrix<double>>>
      d2K; // p x p second derivatives
  // Observation-side linear transform Z and its parameter derivatives
  SparseMatrix<double, 0, int> Z;
  vector<SparseMatrix<double, 0, int>> dZ;
  std::vector<std::vector<SparseMatrix<double, 0, int>>>
      d2Z; // p x p second derivatives
  // Optional: per-parameter fixing mask (K-parameters)
  std::vector<bool> fix_mask_K;

  // Internal solvers for K factorizations
  sparse_llt_solver cholK_solver; // for K or K^T K
  bool llt_inited{false};

  // Global modes (fixed at construction)
  bool analyzed_cholK{false};
  VectorXd trace_vals; // size n_theta_K; tr(K^-1 dK) or NormalEq variant
  bool trace_ready{false};
  MatrixXd HK_trace; // n_theta_K x n_theta_K; H_K trace block
  bool HK_trace_ready{false};

public:
  Operator(const Rcpp::List &operator_list)
      : h(Rcpp::as<VectorXd>(operator_list["h"])),
        n_theta_K(Rcpp::as<int>(operator_list["n_theta_K"])),
        zero_trace(Rcpp::as<bool>(operator_list["zero_trace"])),
        symmetric(Rcpp::as<bool>(operator_list["symmetric"])),
        generic_type(Rcpp::as<string>(operator_list["generic_type"])),
        K(Rcpp::as<SparseMatrix<double>>(operator_list["K"])), dK(n_theta_K) {
    // initial dK
    for (int i = 0; i < n_theta_K; i++) {
      dK[i].resize(h.size(), h.size());
      dK[i].setZero();
    }
    // Initialize Z and dZ
    Z.resize(h.size(), h.size());
    Z.setIdentity();
    dZ.resize(n_theta_K);
    for (int i = 0; i < n_theta_K; ++i) {
      dZ[i].resize(h.size(), h.size());
      dZ[i].setZero();
    }
  }
  virtual ~Operator();

  int get_n_theta_K() const { return n_theta_K; }
  const VectorXd &get_h() const { return h; }
  bool is_symmetric() const { return symmetric; }
  bool is_zero_trace() const { return zero_trace; }

  const SparseMatrix<double> &getK() const { return K; }
  const SparseMatrix<double, 0, int> &get_dK(int i) const { return dK[i]; }
  const SparseMatrix<double, 0, int> &getZ() const { return Z; }
  const SparseMatrix<double, 0, int> &get_dZ(int i) const { return dZ[i]; }
  const SparseMatrix<double, 0, int> &get_d2K(int i, int j) const {
    return d2K[i][j];
  }
  const SparseMatrix<double, 0, int> &get_d2Z(int i, int j) const {
    return d2Z[i][j];
  }
  const MatrixXd &get_HK_trace() const { return HK_trace; }

  // Core builders for K and Z
  // New unified builder: preferred override in subclasses
  virtual void build_KZ(const VectorXd &theta) = 0;

  // Optional analytic derivatives: return true if both dK and/or dZ were set
  virtual bool update_dKdZ(const VectorXd &) {
    return false;
  } // default: no analytic
  // Optional: analytic second derivatives
  virtual bool update_d2Kd2Z(const VectorXd &) { return false; }

  // Optional: per-parameter fixing mask for theta_K (default: none fixed)
  virtual std::vector<bool> get_fix_mask_K() const {
    return std::vector<bool>(n_theta_K, false);
  }

  // Unified updater: update K, Z, (optionally) dK, dZ, factorization, and
  // traces
  virtual void update_all(const VectorXd &theta, const UpdateOptions &opts);

  // Accessors for traces
  const VectorXd &get_trace_trK() const { return trace_vals; }
  bool traces_ready() const { return trace_ready; }
};

class Matern : public Operator {
private:
  SparseMatrix<double, 0, int> G, C, Ci;
  // Basis for log kappa: dense in R (used directly in fractional routines)
  MatrixXd Bk_dense;
  bool stationary{true};
  double alpha;
  VectorXd Cdiag;
  bool fix_alpha{true};
  int m{0}; // rational approximation order (0 = none)
  int dim{2};

public:
  Matern(const Rcpp::List &);

  void build_KZ(const VectorXd &) override;
  int get_alpha() const { return alpha; }
};

// ARMA(p,q) operator: K = G + sum_j phi_j C_j; Z = I + sum_k theta_k L^k
class Arma : public Operator {
private:
  int n;
  int p; // AR order
  int q; // MA order
  // Bases
  std::vector<SparseMatrix<double, 0, int>> Cj; // lag j bases for AR
  SparseMatrix<double, 0, int> G;               // identity
  SparseMatrix<double, 0, int> L; // 1-step lag/shift (subdiagonal 1)
  std::vector<SparseMatrix<double, 0, int>> Lpow; // powers of L
  // fixing masks
  std::vector<bool> fix_phi_mask;   // size p
  std::vector<bool> fix_theta_mask; // size q
  // parameter split: first p are phi, last q are theta
public:
  Arma(const Rcpp::List &);
  void build_KZ(const VectorXd &) override;
  std::vector<bool> get_fix_mask_K() const override {
    std::vector<bool> mask(n_theta_K, false);
    for (int j = 0; j < p; ++j)
      mask[j] = (fix_phi_mask.size() == (size_t)p) ? fix_phi_mask[j] : false;
    for (int k = 0; k < q; ++k)
      mask[p + k] =
          (fix_theta_mask.size() == (size_t)q) ? fix_theta_mask[k] : false;
    return mask;
  }
};

// OU (Ornstein-Uhlenbeck) operator: band matrix with exp(-theta*dt)
// coefficients
class OU : public Operator {
private:
  int n;
  VectorXd dt; // time differences between mesh locations
public:
  OU(const Rcpp::List &);
  void build_KZ(const VectorXd &) override;
};

class Tensor_prod : public Operator {
private:
  std::shared_ptr<Operator> first, second;
  int n_theta_1, n_theta_2;

public:
  Tensor_prod(const Rcpp::List &);

  void build_KZ(const VectorXd &) override;
  bool update_dKdZ(const VectorXd &) override;
};

class Spacetime : public Operator {
private:
  VectorXd Ct_diag, Cs_diag; // not used
  SparseMatrix<double, 0, int> BtCs, Gs, Ct, Cs, Bx, By, S, Bs, Hxx, Hyy, Hxy,
      Hyx;
  // MatrixXd B_gamma_x, B_gamma_y;
  Rcpp::List B_gamma_x_list_input, B_gamma_y_list_input;
  VectorXd theta_gamma_x, theta_gamma_y;
  int n_theta_gamma_x, n_theta_gamma_y;
  double lambda, alpha;
  string method; // galerkin, backward Euler
  bool stabilization, fix_gamma, shared_theta_gamma;
  int nt;
  std::vector<MatrixXd> B_gamma_x_list, B_gamma_y_list;

public:
  Spacetime(const Rcpp::List &);

  void build_KZ(const VectorXd &);
  void update_dK(const VectorXd &);
};

// Bivar
class Bivar : public Operator {
private:
  std::shared_ptr<Operator> first, second;
  int n_theta_1, n_theta_2;
  int n; // dim of K1 and K2 (same)
  bool share_param, fix_theta, use_c_param;
  double bv_theta;

public:
  Bivar(const Rcpp::List &);

  void build_KZ(const VectorXd &) override;

  Matrix2d getD(double, double) const;
  Matrix2d get_dD_theta(double, double) const;
  Matrix2d get_dD_rho(double, double) const;
  Matrix2d get_dD2_theta(double, double) const;
  Matrix2d get_dD2_rho(double, double) const;
};

class Generic : public Operator {
private:
  vector<SparseMatrix<double, 0, int>> matrices;
  vector<string> param_names;
  std::map<string, vector<string>> trans_map;

public:
  Generic(const Rcpp::List &);

  void build_KZ(const VectorXd &theta_K);
  double apply_transform(double value, const string &trans_type) const;
};

class generic_ns : public Operator {
private:
  // Matrices for the model
  vector<SparseMatrix<double, 0, int>> matrices;

  // Position combinations for matrix operations
  vector<vector<int>> position;

  // Parameter mappings
  std::map<string, vector<int>> param_map;

  // Parameter transformations
  std::map<string, string> trans_map;

  // Basis matrices for spatial varying parameters
  std::map<string, MatrixXd> B_theta_K;

public:
  generic_ns(const Rcpp::List &);

  void build_KZ(const VectorXd &theta_K);
  double apply_transform(double value, const string &trans_type) const;
};

// Bivar_normal_ope (theta=0)
class bv_matern : public Operator {
private:
  std::shared_ptr<Matern> first, second;
  int n_theta_1, n_theta_2;
  int n; // dim of K1 and K2 (same)
  bool share_param, fix_theta;
  double dim, alpha1, alpha2, nu1, nu2, bv_theta;

public:
  bv_matern(const Rcpp::List &);

  void build_KZ(const VectorXd &) override;

  Matrix2d getD(double, double) const;
  Matrix2d get_dD_theta(double, double) const;
  Matrix2d get_dD_rho(double, double) const;
  Matrix2d get_dD2_theta(double, double) const;
  Matrix2d get_dD2_rho(double, double) const;
};

// ---- Structure for random effects ----
// U|V ~ N(0, Sigma)
class Randeff : public Operator {
private:
  int n_reff;

public:
  Randeff(const Rcpp::List &);

  void build_KZ(const VectorXd &theta_K) override;
  bool update_dKdZ(const VectorXd &theta_K) override;
};

// for initialize Latent models
class OperatorFactory {
public:
  static std::shared_ptr<Operator> create(const Rcpp::List &operator_in);
};

#endif
