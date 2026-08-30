/*
BlockModel - for ngme each replicate
*/
#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include "include/MatrixAlgebra.h"
#include "include/factor_counters.h"
#include "include/solver.h"
#include "include/timer.h"
#include "latent.h"
#include "model.h"
#include "noise.h"
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using std::vector;

const int BLOCK_FIX_FLAG_SIZE = 7;
enum Block_fix_flag {
  block_fix_beta,
  // noise_parameters
  block_fix_theta_mu,
  block_fix_theta_sigma,
  block_fix_theta_nu,
  block_fix_rho,
  blcok_fix_V,
  // currently not supported
  block_fix_theta_sigma_normal
};

enum precond_type { precond_none, precond_fast, precond_full };

class BlockModel {
protected:
  // W_sizes = row(A1) + ... + row(An)
  // general
  std::mt19937 rng;

  MatrixXd X;
  VectorXd Y;
  int W_sizes, V_sizes; // V_sizes = sum(nrow(K_i))
  string family;

  // Fixed effects and Measurement noise
  VectorXd theta_mu, theta_sigma, theta_nu, rho, beta;

  MatrixXd B_mu, B_sigma, B_nu;
  VectorXd noise_mu, noise_sigma, noise_nu;
  int n_theta_mu, n_theta_sigma, n_theta_nu, n_rho;

  double nu_lower_bound{0.0};

  int n_latent; // how mnay latent model
  int n_obs;    // how many observation

  // n_merr = noise_nu.size + noise_sigma.size + noise_nu.size + rho.size
  // beta.size
  int n_la_params, n_feff, n_merr, n_repl; // number of total replicates(blocks)

  // Correlated measurement error
  bool corr_measure;
  vector<int> cor_rows, cor_cols;
  vector<bool> has_correlation;
  int n_corr_pairs;
  SparseMatrix<double> Q_eps, dQ_eps, sqrt_Rinv;
  // Q_eps = A^T diag(1/SV) A (no correlation)
  // Q_eps = A^T vcov^-1 A (with correlation)
  // vcov = diag(sigma sqrt(V)) R diag(sigma sqrt(V))
  int n_params;

  // fix estimation
  bool fix_flag[BLOCK_FIX_FLAG_SIZE]{0};
  std::vector<bool> fix_theta_sigma_vec;

  // controls
  int n_gibbs;
  bool debug, reduce_var;
  bool robust{false};
  double reduce_power, threshold;

  SparseMatrix<double> A, K, Q, QQ, pmat, pmat_inv;

  vector<std::shared_ptr<Latent>> latents;
  VectorXd p_vec, a_vec, b_vec, noise_V, noise_prevV;
  // double nu {1};

  // optimize related
  VectorXd stepsizes;
  int counting{0};
  VectorXd indicate_threshold, steps_to_threshold;
  int curr_iter; // how many times set is called.

  bool all_gaussian, rao_blackwell,
      use_iterative_solver; // No need for gibbs sampling
  std::vector<std::string> par_names;
  VectorXd rb_trace_noise_sigma;
  // Per-parameter variance contributed by the Hutchinson trace estimators
  // (probe variance / N_iter), in the same layout as the per-latent traces.
  std::vector<VectorXd> rb_probe_var_K_latent, rb_probe_var_sigma_latent;
  VectorXd rb_probe_var_noise_sigma;
  // Running (Welford) variance of the gradient across iterations, used to size
  // the probe budget against the Gibbs noise rather than fixing it a priori.
  // Lag-1 differences rather than a plain running variance: during optimization
  // the gradient's mean is moving, and a raw variance charges that drift to the
  // sampling noise, which makes the Gibbs term look larger than it is and
  // under-provisions probes. E[(g_t - g_{t-1})^2] / 2 is unaffected by a slowly
  // moving mean.
  VectorXd grad_prev_, grad_diff_sq_;
  VectorXd probe_var_ewma_;  // smoothed raw probe variance, parameter layout
  long grad_run_n_{0};
  int trace_adapt_countdown_{0};
  int last_trace_N_{0};
  int selinv_state_{-1}; // -1 undecided, 0 probing, 1 selected inverse
  double selinv_max_fill{4.0};
  bool trace_adapt{false};
  double trace_adapt_frac{0.1};
  int trace_adapt_every{50}, trace_adapt_min{5}, trace_adapt_max{200};
  // Store per-latent RB trace terms at Block level
  std::vector<Eigen::VectorXd> rb_trace_K_latent;
  std::vector<Eigen::VectorXd> rb_trace_sigma_latent;
  // Option: use conditional mean of W instead of sampling when estimating
  bool use_cond_W{false};

  // Gradient covariance matrix storage
  MatrixXd grad_covariance;

  // Cached preconditioners (averages from the most recent grad loop)
  MatrixXd last_precond; // cached preconditioner for last strategy
  bool last_precond_valid{false};

  // Cached gradient from the last compute pass
  Eigen::VectorXd last_gradient;
  bool last_grad_valid{false};

  // priors
  std::vector<string> prior_beta_type;
  std::vector<VectorXd> prior_beta_param;
  std::vector<string> prior_beta_target;
  string prior_mu_type, prior_sigma_type, prior_nu_type;
  VectorXd prior_mu_param, prior_sigma_param, prior_nu_param;
  string prior_mu_target{"coef"}, prior_sigma_target{"coef"},
      prior_nu_target{"coef"};

  // For computing RB gradient_K
  vector<vector<SparseMatrix<double>>> block_dK;

  sparse_llt_solver chol_QQ;

  // caching of QQ, its Cholesky factor, and the A*Z block matrix
  // QQ = K' diag(1/SV) K + (AZ)' D (AZ) only depends on the current
  // parameters (through K, Z, sigma, noise_sigma, Q_eps) and on the latent /
  // measurement variance vectors V. None of those change between the Gibbs
  // draws of a Gaussian model, so assembling and refactorizing QQ once per
  // optimizer iteration is enough. QQ_valid records whether the cached QQ and
  // chol_QQ still match the current state; every mutator that can invalidate
  // them calls invalidate_QQ() / invalidate_AZ().
  bool QQ_valid{false};
  // Sparsity pattern QQ had when chol_QQ.analyze() was last run. The symbolic
  // phase only has to be redone when the pattern actually changes, which for
  // rational (fractional) approximations happens when the smoothness crosses
  // an integer and the operator gains or loses factors.
  bool QQ_analyzed{false};
  std::vector<int> QQ_pat_outer, QQ_pat_inner;
  // Remember / compare the sparsity pattern chol_QQ was last analyzed for.
  void record_QQ_pattern();
  bool QQ_pattern_changed() const;
  // AZ = [A_1 Z_1, ..., A_L Z_L] is rebuilt three times per Gibbs draw in the
  // uncached code; it only depends on the parameters, so cache it alongside QQ.
  mutable SparseMatrix<double> AZ_cached;
  mutable bool AZ_valid{false};
  // Measurement-side pieces. H = sqrt(D) A Z (times sqrt_Rinv in the
  // correlated case) and the measurement block of QQ -- H'H, or
  // (AZ)' Q_eps (AZ). This depend on AZ and on the measurement precision only,
  // never on the latent V. They are therefore constant across the Gibbs draws
  // of any model with gaussian measurement noise, even when the latent noise
  // is non-gaussian, where the latent block K' diag(1/SV) K does change every
  // draw. Splitting the two halves lets the expensive n_obs-sized product be
  // reused for NIG/GAL/t latent fields.
  mutable SparseMatrix<double> sqrt_AtSVA_cached;
  mutable bool sqrt_AtSVA_valid{false};
  mutable SparseMatrix<double> QQ_measure;
  mutable bool QQ_measure_valid{false};
  // G = sqrt(1/SV) K, the latent half of the rMVN draw. Same dependencies as
  // the latent block of QQ (K and the latent V), so it is rebuilt on exactly
  // the same occasions
  mutable SparseMatrix<double> G_cached;
  mutable bool G_valid{false};

  // clock
  std::chrono::milliseconds sampling_time{0}, update_time{0};

public:
  // BlockModel() {}
  BlockModel(const Rcpp::List &block_model, unsigned long seed);
  ~BlockModel() = default;

  /* Gibbs Sampler */
  void burn_in(int);

  int get_n_obs() const { return n_obs; }
  void sampleW_VY(bool burn_in = false);

  bool is_all_gaussian() const { return all_gaussian; }

  void sample_cond_V() {
    if (n_latent > 0) {
      for (int i = 0; i < n_latent; i++) {
        if (latents[i]->V_may_change())
          invalidate_QQ();
        latents[i]->sample_cond_V();
      }
    }
  }

  void sample_uncond_V() {
    if (n_latent > 0) {
      for (int i = 0; i < n_latent; i++) {
        if (latents[i]->V_may_change())
          invalidate_QQ();
        latents[i]->sample_uncond_V();
      }
    }
  }

  // Mark the latent-side caches out of date. Called from the latent V
  // samplers, which change SV and hence the latent block of QQ and G.
  void invalidate_QQ() {
    QQ_valid = false;
    G_valid = false;
  }
  // The measurement precision moved (noise_sigma, noise_V, Q_eps or rho), so
  // H and the measurement block of QQ have to be rebuilt as well.
  void invalidate_measurement() {
    sqrt_AtSVA_valid = false;
    QQ_measure_valid = false;
    QQ_valid = false;
  }
  // The parameters moved, so K, Z and everything downstream is stale.
  void invalidate_AZ() {
    AZ_valid = false;
    G_valid = false;
    invalidate_measurement();
  }

  // AZ = [A_1 Z_1, ..., A_L Z_L], rebuilt only when the parameters changed.
  const SparseMatrix<double> &get_AZ() const;
  // Measurement block of QQ, rebuilt only when the measurement precision or
  // the parameters changed.
  const SparseMatrix<double> &get_QQ_measure() const;
  // G = sqrt(1/SV) K, with the same un-clamped 1/SV that sampleW_VY() uses.
  const SparseMatrix<double> &get_G(const VectorXd &inv_SV) const;

  // Reassemble + refactorize QQ only if the cache is stale.
  void ensure_QQ() {
    if (!QQ_valid || ngme_counters::cache_disabled())
      update_QQ();
  }

  void update_QQ();
  void update_Q_eps(double rho);
  const SparseMatrix<double> &get_sqrt_AtSVA() const;

  void setW(const VectorXd &);
  void setPrevW(const VectorXd &);
  void setPrevV(const VectorXd &);
  void set_cond_W(const VectorXd &);

  /* Optimizer related */
  int get_n_params() const { return n_params; }

  VectorXd get_parameter();
  void set_parameter_and_update(const VectorXd &, bool with_precond);

  // Accessors after compute
  MatrixXd get_preconditioner();
  const VectorXd &get_gradient() const {
    if (!last_grad_valid) {
      throw std::runtime_error("last_gradient is not valid");
    }
    return last_gradient;
  }
  MatrixXd get_grad_covariance() const { return grad_covariance; }

  // Unified compute for both gradient and (optionally) preconditioner
  void compute_grad_and_hessian(bool with_precond, double eps);

  /* Aseemble */
  void assemble() {
    int nrow = 0;
    int ncol = 0;
    for (vector<std::shared_ptr<Latent>>::iterator it = latents.begin();
         it != latents.end(); it++) {
      setSparseBlock(&K, nrow, ncol, (*it)->getK());
      // setSparseBlock(&dK,  n, n, (*it)->get_dK());
      // setSparseBlock(&d2K, n, n, (*it)->get_d2K());
      nrow += (*it)->get_V_size();
      ncol += (*it)->get_W_size();
    }
  }

  // assemble dK, and dK_sigma
  void assemble_dK();

  // tr(QQ^-1 dK^T diag(1/SV) K)
  void compute_rb_trace();
  void adapt_trace_probes();
  // tr(QQ^{-1} T): exact via the selected inverse when the factor is cheap
  // enough, Hutchinson otherwise. Records the probe variance either way (zero
  // for the exact path) so the probe adapter keeps working.
  double qq_trace(const Eigen::SparseMatrix<double, 0, int> &T, double &probe_var);
  int get_trace_N() const { return chol_QQ.get_N_iter(); }

  // return mean = mu*(V-h)
  VectorXd getMean() const {
    VectorXd mean(V_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_V_size();
      mean.segment(pos, size) = (*it)->getMean();
      pos += size;
    }
    return mean;
  }

  VectorXd getV() const {
    VectorXd V(V_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_V_size();
      V.segment(pos, size) = (*it)->getV();
      pos += size;
    }

    return V;
  }

  VectorXd getPrevV() const {
    VectorXd V(V_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_V_size();
      V.segment(pos, size) = (*it)->getPrevV();
      pos += size;
    }

    return V;
  }

  // return sigma * V
  VectorXd getSV() const {
    VectorXd SV(V_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_V_size();
      SV.segment(pos, size) = (*it)->getSV();
      pos += size;
    }

    return SV;
  }

  VectorXd getW() const {
    VectorXd W(W_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_W_size();
      W.segment(pos, size) = (*it)->getW();
      pos += size;
    }
    return W;
  }

  VectorXd get_cond_W() const {
    VectorXd W(W_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_W_size();
      W.segment(pos, size) = (*it)->get_cond_W();
      pos += size;
    }
    return W;
  }

  VectorXd getPrevW() const {
    VectorXd W(W_sizes);
    int pos = 0;
    for (vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin();
         it != latents.end(); it++) {
      int size = (*it)->get_W_size();
      W.segment(pos, size) = (*it)->getPrevW();
      pos += size;
    }
    return W;
  }

  VectorXd get_residual(bool rao_blackwell = false) const {
    // Compute sum_i A_i * Z_i * W_i (or cond_W_i) across latents
    if (n_latent > 0) {
      VectorXd AZW = VectorXd::Zero(n_obs);
      for (int i = 0; i < n_latent; ++i) {
        const auto &Ai = latents[i]->getA();
        const auto &Zi = latents[i]->getZ();
        VectorXd Wi =
            rao_blackwell ? latents[i]->get_cond_W() : latents[i]->getW();
        AZW.noalias() += Ai * (Zi * Wi);
      }
      return Y - AZW - X * beta -
             (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
    } else {
      return Y - X * beta -
             (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
    }
  }

  // residual_part = Y - X beta - (1 - V) mu
  VectorXd get_residual_part() const {
    return Y - X * beta -
           (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
  }

  void sample_cond_noise_V(bool posterior = true);

  // for updating hessian
  vector<VectorXd> get_VW() const {
    vector<VectorXd> ret(3);
    ret[0] = noise_V;
    ret[1] = getV();
    ret[2] = getW();
    return ret;
  }

  void set_prev_VW(const vector<VectorXd> &VW) {
    noise_prevV = VW[0];
    setPrevV(VW[1]);
    setPrevW(VW[2]);
  }

  // Methods for managing Gibbs samples for preconditioner
  // --------- Fixed effects and Measurement error  ------------
  VectorXd grad_beta();

  VectorXd get_theta_merr() const;
  VectorXd grad_theta_mu();
  VectorXd grad_theta_sigma();
  VectorXd grad_theta_merr();
  void set_theta_merr(const VectorXd &theta_merr);

  // get length of W,V of iterations
  Rcpp::List sampling(int n, int n_burnin, bool posterior,
                      const SparseMatrix<double> &A);

  Rcpp::List sampling(int n, int n_burnin, bool posterior) {
    return sampling(n, n_burnin, posterior, A);
  }

  Rcpp::List output() const;
  std::vector<std::string> get_par_names() const { return par_names; }

  static double th2rho(double th) { return (-1 + 2 * exp(th) / (1 + exp(th))); }
  static double rho2th(double r) { return (log((-1 - r) / (-1 + r))); }

  // drho / dtheta
  // double dtheta_rho(double th) const {return 2 * exp(th) / pow(1+exp(th),
  // 2);}
  static double dtheta_th(double rho) { return (1 - rho * rho) / 2; }

  // numerical Hessian path removed; all measurement/fixed-effects Hessians are
  // analytic now

  double log_likelihood() { return 0.0; }
};

#endif
