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
  // tr(K^{-1} dK). Cleared by callers that do not need it (theta_K fixed, or a
  // zero-trace operator), which lets update_all skip factorizing K altogether.
  bool compute_trace{true};
  // Seed for the Hutchinson probe vectors.
  unsigned int trace_seed{0};
  // Deprecated / no longer read: the symbolic analyze() is now re-run exactly
  // when K's sparsity pattern changes, which is strictly better than an
  // unconditional refresh (analyze() is a pure function of the pattern).
  bool robust_reanalyze{false};
  bool prefer_analytic_dK{true};
  bool prefer_analytic_dZ{true};
  // true so an operator that HAS a closed form for d2K uses it. Every operator
  // without one returns false from update_d2Kd2Z() and still gets the numeric
  // route, so this changes nothing for them.
  bool prefer_analytic_d2K{true};
  bool prefer_analytic_d2Z{false};
  int n_trace_iter{8};
  // Fill above which the exact selected inverse is not considered at all. Only
  // a rail here, not the decision: control_opt's selinv_max_fill_k feeds this
  // and is Inf by default, so the count below normally decides. Distinct from
  // BlockModel's selinv_max_fill, which IS the decision for the block
  // precision.
  double selinv_max_fill{4.0};
  // How much dearer one exact selected inverse may be than the probes it
  // replaces and still be preferred; the exact route carries no estimation
  // variance, so it is worth somewhat more than its bare cost.
  double selinv_cost_ratio{2.0};
  // Echoes the latent's debug flag so the route chosen can be printed.
  bool debug{false};
  int solver_type{0};
  // 0 = LU of K, 1 = Cholesky of K^T K (non-symmetric operators only)
  int nonsym_solver{0};
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
  // Sparsity pattern K had when cholK_solver.analyze() was last run. For most
  // operators the pattern is fixed and the symbolic phase runs exactly once;
  std::vector<int> K_pat_outer, K_pat_inner;
  void record_K_pattern();
  bool K_pattern_changed() const;
  // eps the current update_all() would use for numeric differencing, so an
  // override of update_dKdZ can reproduce the base class's difference exactly
  // instead of inventing its own step.
  double last_eps_dK_{1e-4};
  // The options of the update_all() call currently in progress, so that
  // update_dKdZ / update_d2Kd2Z can see what the caller actually asked for.
  // Valid only inside update_all().
  const UpdateOptions *cur_opts_{nullptr};
  VectorXd trace_vals; // size n_theta_K; tr(K^-1 dK) or NormalEq variant
  bool trace_ready{false};
  MatrixXd HK_trace; // n_theta_K x n_theta_K; H_K trace block
  bool HK_trace_ready{false};
  // -1 undecided, 0 probes, 1 exact selected inverse. Decided once: it turns on
  // the fill of K's factor, and K's pattern does not move during a fit.
  int selinv_state_{-1};

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
  // d2K / d2Z are allocated only when second derivatives are actually asked
  // for, so these have to answer for an operator that was never asked. Callers
  // test the returned matrix (rows() > 0) to decide whether a second-derivative
  // term exists; returning an empty one keeps that test correct rather than
  // indexing an empty vector.
  const SparseMatrix<double, 0, int> &get_d2K(int i, int j) const {
    static const SparseMatrix<double, 0, int> empty;
    if (i < 0 || j < 0 || (int)d2K.size() <= i || (int)d2K[i].size() <= j)
      return empty;
    return d2K[i][j];
  }
  const SparseMatrix<double, 0, int> &get_d2Z(int i, int j) const {
    static const SparseMatrix<double, 0, int> empty;
    if (i < 0 || j < 0 || (int)d2Z.size() <= i || (int)d2Z[i].size() <= j)
      return empty;
    return d2Z[i][j];
  }
  const MatrixXd &get_HK_trace() const { return HK_trace; }

  // Core builders for K and Z
  // New unified builder: preferred override in subclasses
  virtual void build_KZ(const VectorXd &theta) = 0;

  // Structural shortcut for the traces. An operator whose K has exploitable
  // structure can fill trace_vals (and HK_trace, when opts.compute_HK_trace)
  // itself and return true; update_all then skips the factorization of K
  // entirely, since cholK_solver exists for nothing else. Default: no
  // shortcut, use the generic Hutchinson path.
  virtual bool compute_traces_structured(const VectorXd & /*theta*/,
                                         const UpdateOptions & /*opts*/) {
    return false;
  }

protected:
  // Exact O(n) traces for a triangular K, see operator.cpp. Returns false if K
  // is not triangular (or has a zero on the diagonal), leaving the caller to
  // fall back on the factorization.
  bool try_triangular_traces(const UpdateOptions &opts, bool want_trace,
                             bool want_HK);

public:

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
  // Closed-form derivatives for the integer cases. kappa = exp(B_kappa theta)
  // enters K only through a diagonal, so every derivative is one diagonal
  // rescaling plus at most two sparse products.
  bool update_dKdZ(const VectorXd &) override;
  bool update_d2Kd2Z(const VectorXd &) override;

private:
  // A = G + (C diag(kappa^2)  or  Dk C Dk), the alpha = 2 operator and the
  // inner factor of the alpha = 4 one. Also returns the per-parameter
  // derivatives of A, which both orders need.
  bool matern_dA(const VectorXd &theta_K, SparseMatrix<double> &A,
                 std::vector<SparseMatrix<double>> &dA,
                 std::vector<std::vector<SparseMatrix<double>>> *d2A) const;
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
  // True once the factors have been brought up to date at the current theta
  // during THIS update_all, so the three places that need their derivatives
  // and traces share one factor update instead of repeating it. Cleared in
  // build_KZ(), which runs at the start of every update_all.
  bool factors_current_{false};
  // Build the per-factor options and run both factors' update_all.
  bool update_factors(const VectorXd &theta, const UpdateOptions &opts);

public:
  Tensor_prod(const Rcpp::List &);

  void build_KZ(const VectorXd &) override;
  // Analytic derivatives from the Kronecker structure:
  //     dK/dtheta_1j = (dK_1/dtheta_1j) (x) K_2
  //     dK/dtheta_2j = K_1 (x) (dK_2/dtheta_2j)
  // Without this the base class falls back to numeric differencing, which
  // rebuilds the whole Kronecker product once per parameter per iteration.
  bool update_dKdZ(const VectorXd &) override;
  // Second derivatives from the same Kronecker structure:
  //     d2K/dtheta_1j dtheta_1k = (d2K_1/dtheta_1j dtheta_1k) (x) K_2
  //     d2K/dtheta_1j dtheta_2k = (dK_1/dtheta_1j) (x) (dK_2/dtheta_2k)
  //     d2K/dtheta_2j dtheta_2k = K_1 (x) (d2K_2/dtheta_2j dtheta_2k)
  // Without this the base class differences the WHOLE Kronecker product four
  // times per parameter pair.
  bool update_d2Kd2Z(const VectorXd &) override;
  bool compute_traces_structured(const VectorXd &,
                                 const UpdateOptions &) override;
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
  // True when every entry of B_gamma_x_list is the same matrix, and likewise
  // for y. The advection field is then the same at every time node, so all
  // nt - 1 spatial blocks of K are the same matrix and build_KZ() assembles one
  // instead of nt - 1. Decided once in the constructor: the B lists are fixed
  // data, not parameters, so this cannot change as theta moves.
  bool gamma_time_invariant{false};
  // BtCs with its first block row removed. That row is the rw1 operator's
  // trapezoid row, which couples time slice 1 to every other slice; dropping it
  // is what leaves K block lower-bidiagonal and lets slice 1 carry its own
  // stationary block instead.
  SparseMatrix<double, 0, int> BtCs_st;
  int ns_{0};                  // spatial dimension = Cs.rows()
  bool stationary_init{true};  // give slice 1 the stationary distribution
  // Put the while cc factor on the temporal term instead of splitting it
  // between temporal and spatial, to reduce confounding with sigma.
  bool cc_variance_free{false};
  bool block_trace_checked_{false};

  // Selected-inverse route for the diagonal-block traces. tr(M^-1 B) needs
  // M^-1 only where B is nonzero, so one selected inverse per distinct block
  // serves every parameter . One solver per distinct block, kept across
  // iterations so the symbolic phase is paid once: only the values move.
  std::vector<std::unique_ptr<sparse_llt_solver>> blk_solver_;
  // nnz the block had when that solver was last analyzed; -1 = never.
  std::vector<long long> blk_nnz_;
  // -1 undecided, 0 ruled out (non-symmetric block, failed factorization, or a
  // derivative entry outside the factor pattern) -- once ruled out, stay on the
  // LU path rather than retrying the same failure every iteration.
  int blk_selinv_state_{-1};
  // Debug counter for NGME_SPACETIME_TRACE_CHECK only.
  // LU per distinct diagonal block, built only when the Hessian traces need
  // dense solves and the selected-inverse path (which has no LU) is in use.
  // By pointer: Eigen's SparseLU is neither copyable nor movable, so a plain
  // vector of them cannot be resized.
  std::vector<std::unique_ptr<Eigen::SparseLU<SparseMatrix<double>>>> hk_lu_;
  int selinv_calls_{0};
  int selinv_fail_{0};

  // tr(M_b^-1 B) via the selected inverse of block b. False when the route is
  // unavailable, in which case the caller must fall back.
  bool block_trace_selinv(int b, const SparseMatrix<double, 0, int> &B,
                          double &out);

  // ---- closed-form derivative helpers, see tensorprod.cpp ----
  // out = base_coef * (BtCs_st or BtCs) + blockdiag(0, blk...) + pad(first).
  // Every derivative of K has this shape, because K itself does.
  void assemble_shaped(SparseMatrix<double> &out, double base_coef,
                       const std::vector<SparseMatrix<double>> &blk,
                       bool uniform, const SparseMatrix<double> &first,
                       bool has_first) const;
  // dLs/dtheta_m and d2Ls/dtheta_m dtheta_n for the interior block at time node
  // i. theta_0 (cc) does not enter Ls, so m, n >= 1. False means identically
  // zero, which the caller can then skip rather than assemble.
  bool dLs_block(int m, int i, double k2, SparseMatrix<double> &out) const;
  bool d2Ls_block(int m, int n, int i, double k2,
                  SparseMatrix<double> &out) const;
  // Everything both derivative routines need from theta, gathered once. False
  // for a configuration the closed form does not cover.
  bool analytic_state(const VectorXd &theta_K, double &c, double &k2, double &a,
                      double &s, double &q, bool &uniform,
                      std::vector<SparseMatrix<double>> &Ls, bool &has_stat,
                      SparseMatrix<double> &M, double &mm, double &nn,
                      double &u) const;

public:
  Spacetime(const Rcpp::List &);

  void build_KZ(const VectorXd &);
  void update_dK(const VectorXd &);
  // Exact traces from the diagonal blocks alone. K is block lower-bidiagonal
  // once the stationary initial condition is in place, so this needs neither
  // probes nor a factorization of the full (nt*ns) x (nt*ns) operator.
  bool compute_traces_structured(const VectorXd &,
                                 const UpdateOptions &) override;

  // Closed-form derivatives. Numeric differencing calls build_KZ() once per
  // parameter for dK and FOUR times per parameter PAIR for d2K, and every one
  // of those assembles the full (nt*ns) operator. Both return false for a
  // configuration they do not cover, leaving the numeric route in place.
  bool update_dKdZ(const VectorXd &) override;
  bool update_d2Kd2Z(const VectorXd &) override;
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
  // K is a linear combination of fixed matrices whose coefficients are
  // products of scalar transforms of theta, so every derivative is available
  // in closed form.
  bool update_dKdZ(const VectorXd &) override;
  bool update_d2Kd2Z(const VectorXd &) override;

private:
  // Per-matrix coefficient factors for each parameter, and their derivatives:
  // f[p][i] = T(theta_p, trans[p][i]), with f1 and f2 the 1st and 2nd
  // derivatives. A parameter absent from trans_map contributes a constant 1.
  void coef_factors(const VectorXd &theta_K, std::vector<VectorXd> &f,
                    std::vector<VectorXd> &f1, std::vector<VectorXd> &f2) const;
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

// ---- VAR(1) bivariate operator with Cayley reparameterization ----
// theta_K = (p1, p2, p3, p4): four unconstrained parameters.
// A = Cayley(J, L) guarantees spectral radius < 1.
// K = M0 + a11*M11 + a22*M22 + a12*M12 + a21*M21.
// Numerical dK is computed automatically by update_all(). Thread-safe.
class RCallback : public Operator {
private:
  using SM = SparseMatrix<double, 0, int>;
  SM M0, M11, M22, M12, M21;
  double eps_cayley;

public:
  RCallback(const Rcpp::List &);
  void build_KZ(const VectorXd &theta_K) override;
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
