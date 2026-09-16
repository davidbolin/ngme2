/*
BlockModel - for ngme each replicate
*/
#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include "include/MatrixAlgebra.h"
#include "include/factor_counters.h"
#include "include/nig_std.h"
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
  int nig_param_std{0};
  double reduce_power, threshold;

  SparseMatrix<double> A, K, Q, QQ, pmat, pmat_inv;

  // Per-observation weight on the MEASUREMENT precision, all ones normally.
  // Zeroing an entry removes that observation from QQ and from M exactly, which
  // is how the exact leave-group-out chains drop a fold without rebuilding the
  // model: the latent block, the operator and the symbolic factorisation are
  // untouched, and dropping rows can only shrink QQ's pattern.
  VectorXd obs_weight;
  // Saved initial state; see snapshot_state().
  std::vector<VectorXd> snap_W, snap_V;
  VectorXd snap_noise_V;
  bool snap_taken{false};

  vector<std::shared_ptr<Latent>> latents;
  VectorXd p_vec, a_vec, b_vec, noise_V, noise_prevV;
  // double nu {1};

  // optimize related
  VectorXd stepsizes;
  int counting{0};
  VectorXd indicate_threshold, steps_to_threshold;
  int curr_iter; // how many times set is called.

  bool all_gaussian, rao_blackwell;
  std::vector<std::string> par_names;
  VectorXd rb_trace_noise_sigma;
  // Per-parameter variance contributed by the Hutchinson trace estimators
  // (probe variance / N_iter), in the same layout as the per-latent traces.
  std::vector<VectorXd> rb_probe_var_K_latent, rb_probe_var_sigma_latent;
  VectorXd rb_probe_var_noise_sigma;
  // Running (Welford) variance of the gradient across iterations, used to size
  // the probe budget against the Gibbs noise rather than fixing it a priori.
  VectorXd grad_prev_, grad_diff_sq_;
  VectorXd probe_var_ewma_;  // smoothed raw probe variance, parameter layout
  long grad_run_n_{0};
  // The budget this block would like next, or -1 if it has nothing to say.
  // adapt_trace_probes() writes it; it is NOT applied here. Parallel chains
  // must all probe at the SAME budget: R_hat compares chains against each
  // other, and two chains estimating the same gradient at different probe
  // counts have different estimator variances, so their spread no longer
  // measures what the diagnostic assumes it measures. Chains that adapted
  // independently also do unequal work behind a statically scheduled barrier.
  // The driver collects these suggestions where the chains have just joined --
  // the convergence checkpoint -- and hands one budget back to all of them.
  int suggested_trace_N_{-1};
  int selinv_state_{-1}; // -1 undecided, 0 probing, 1 selected inverse
  double selinv_max_fill{4.0};
  // Probe budget for the expected-information estimate used by the
  // preconditioner. It was taking the trace probe budget, which is sized
  // against gradient variance -- a different target -- so it has its own.
  // NA on the R side means follow n_trace_iter, as before.
  int n_fisher_probes_{-1};
  bool trace_adapt{false};
  double trace_adapt_frac{0.1};
  // "share" = size the budget so the trace estimator carries trace_adapt_frac
  // of the gradient variance (the original rule). "cost" = size it to minimise
  // the variance of the Polyak-Ruppert average per unit of work, and only from
  // the polish phase on. See suggest_trace_N_cost().
  bool trace_adapt_cost_rule{false};
  // Probe-proportional share of one gradient computation, measured once when
  // the polish begins and then held. cost_ratio_ is a/b: how many probe
  // columns cost as much as everything else in a pass. -1 = not yet measured.
  double cost_ratio_{-1.0};
  // Fill of QQ's Cholesky factor, as the selinv gate measured it once for its
  // own decision. Reused as the cost ratio; see begin_polish_trace_rule().
  double qq_fill_{-1.0};
  // True once the polish has begun. It changes what the probe block has to
  // deliver: the search under the cost rule does not adapt, so it does not need
  // the probe variance and can run one sign draw over the colouring -- the
  // cheapest form of probing there is. The polish does adapt, and a spread can
  // only be measured across whole replicates, so it needs two.
  bool in_polish_{false};
  // Structure the Hutchinson probes against a colouring of the graph of QQ
  // rather than drawing them densely; see include/probing.h. Off restores the
  // dense Rademacher probes exactly. trace_probing_max_dist caps the colouring
  // distance the solver may consider -- the probe count is capped by the budget
  // either way, so this only bounds the work of looking for a colouring.
  bool trace_probing{true};
  int trace_probing_max_dist{4};
  // The most the probe budget may be multiplied by in order to reach the
  // smallest budget at which probing engages at all. 1 leaves the budget alone.
  // A colouring costs one probe per colour, and that count is a property of the
  // mesh, not of the budget, so on a well-connected graph the smallest usable
  // budget can sit above a small default -- and probing then never engages,
  // however much better it would spend the money.
  double trace_probing_raise_budget{1.0};
  // What the last [probing] debug line said, so the next one is printed only
  // when the configuration actually moved. The budget adapts during a fit and
  // the colouring distance moves with it, so a single line at the start would
  // describe something that is no longer true.
  int probing_reported_dist_{-1}, probing_reported_reps_{-1};
  // How many times the colouring has been re-sourced. See setup_qq_probing().
  int qq_probing_setups_{0};
  // The cadence of budget updates belongs to the driver, which applies one
  // agreed budget to every chain at the convergence checkpoints.
  int trace_adapt_min{5}, trace_adapt_max{200};
  // The operator-side probe block is resized by set_N_iter, and the trace code
  // caches a factored probe block across Gibbs passes. Changing the budget in
  // the middle of a gradient computation therefore leaves that cache sized for
  // the old budget. The new value is parked here and applied at the top of the
  // next computation, where nothing is in flight.
  int pending_k_budget_{-1};
  // Whether trace_adapt also moves the operator-side probe budget. The two were
  // never in step: the QQ budget is adapted while the operator's is fixed at its
  // construction value, although both feed the gradient. Off by default, so the
  // adaptation behaves as before unless asked otherwise.
  bool trace_adapt_k{false};
  // Store per-latent RB trace terms at Block level
  std::vector<Eigen::VectorXd> rb_trace_K_latent;
  std::vector<Eigen::VectorXd> rb_trace_sigma_latent;
  // Option: use conditional mean of W instead of sampling when estimating
  bool use_cond_W{false};

  // Gradient covariance matrix storage
  MatrixXd grad_covariance;
  // Preconditioner for the measurement theta_sigma, control_opt(precond_meas_sigma):
  // 0 = auto (Fisher for non-Gaussian measurement noise), 1 = always Fisher,
  // 2 = complete-data Hessian (non-Gaussian noise falls back to Fisher).
  int precond_meas_sigma_{0};
  // Refresh the cached Fisher block every this many iterations, or sooner once
  // theta_sigma has moved by more than 0.05.
  int fisher_refresh_every_{10};
  // The cached Fisher block, with the iteration and theta_sigma it was built at.
  MatrixXd fisher_sigma_cache_;
  int fisher_cache_iter_{-1000000};
  VectorXd fisher_cache_theta_;

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
  // QQ = Q + QQ_measure is a sparse-sparse ADD, and Eigen recomputes the union
  // pattern and reallocates on every call even though both operands' patterns
  // are fixed within a fit. The nnz each operand had when QQ's pattern was last
  // established; while they hold, the add is a value-only merge into QQ's
  // existing storage.
  long long qq_map_q_nnz_{-1}, qq_map_meas_nnz_{-1};
  void record_qq_add_pattern(const SparseMatrix<double> &Qm,
                             const SparseMatrix<double> &Me);
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

  // The probe budget this block would like next, or -1 if it has none to offer.
  int get_suggested_trace_N() const { return suggested_trace_N_; }
  // Measure the probe/iteration cost ratio and switch the budget rule over to
  // the cost form. Called once, where the polish begins. Returns false if the
  // measurement is not usable, leaving the budget where it is.
  bool begin_polish_trace_rule();
  // Adopt a budget chosen for every chain at once. Applying it is separated
  // from suggesting it so the value can be agreed across chains first; see the
  // note on suggested_trace_N_.
  void apply_trace_N(int N);

  // The drawn W is NOT redundant, even for an all-Gaussian model under
  // Rao-Blackwellisation. The gradient takes cond_W, but the preconditioner's
  // Hessian takes the DRAW, and it must: H_K is QUADRATIC in W, so
  // E[H(W)|Y] != H(E[W|Y]). The difference is tr(A Cov(W|Y)), which the
  // gradient handles with its own RB trace terms and the Hessian has none of.
  void sampleW_VY(bool burn_in = false);
  // Draw W from its prior given V (no conditioning on Y).
  void sampleW_V();

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
  // Point chol_QQ's probes at the graph of QQ (or take them off it when
  // trace_probing is off). Called wherever the symbolic phase runs, since that
  // is exactly where the pattern the colouring describes can have moved.
  void setup_qq_probing();
  // Print the probe structure under debug, but only when it has changed.
  void report_probing();
  // Budget from the cost rule; -1 when it cannot be formed yet.
  int suggest_trace_N_cost();
  // Refresh the smoothed per-parameter probe variances both rules read.
  void update_probe_var_stats();
  double qq_trace(const Eigen::SparseMatrix<double, 0, int> &T, double &probe_var);
  // tr(QQ^-1 A^T diag(d) B). Same quantity as qq_trace() on the assembled
  // triple product, but it hands the three factors to the probe estimator
  // instead of multiplying them out first -- the product is as dense as QQ and
  // is built once per parameter per iteration. Falls back to assembling it when
  // the selected inverse is in use, since that route reads T's entries.
  // Variant taking a pre-computed B*QU, so a loop over parameters that all
  // share the same B pays for that product once. Falls back to the ordinary
  // path whenever the selected inverse is in play (which needs the assembled
  // product anyway).
  double qq_trace_factored_shared(const Eigen::SparseMatrix<double, 0, int> &A,
                                  const Eigen::VectorXd &d,
                                  const Eigen::SparseMatrix<double, 0, int> &B,
                                  const Eigen::MatrixXd &BQU, bool have_BQU,
                                  double &probe_var);
  double qq_trace_factored(const Eigen::SparseMatrix<double, 0, int> &A,
                           const Eigen::VectorXd &d,
                           const Eigen::SparseMatrix<double, 0, int> &B,
                           double &probe_var);
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
  // noise_sigma^-2 / noise_V, with dropped observations zeroed. Every use of
  // the measurement precision during sampling goes through here.
  VectorXd meas_prec() const {
    VectorXd p =
        noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
    if (obs_weight.size() == p.size())
      p = p.cwiseProduct(obs_weight);
    return p;
  }

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
  // Expected information of the marginal likelihood (given V) for the free
  // measurement theta_sigma; uncorrelated measurement noise only.
  MatrixXd fisher_theta_sigma();
  // Mode of the standardised NIG coordinates (nig_std.h) the optimiser uses
  // for the measurement noise; 0 when they do not apply.
  int merr_nig_mode() const;
  VectorXd grad_theta_merr();
  void set_theta_merr(const VectorXd &theta_merr);

  // get length of W,V of iterations
  Rcpp::List sampling(int n, int n_burnin, bool posterior,
                      const SparseMatrix<double> &A);

  Rcpp::List sampling(int n, int n_burnin, bool posterior) {
    return sampling(n, n_burnin, posterior, A);
  }

  // Leave-group-out cross-validation along a SINGLE full-data Gibbs chain.
  //
  // For every retained draw and every group I this records the mean and
  // covariance of eta_I = (AZ W)_I under p(W | y_-I, V), obtained by a
  // rank-|I| downdate of the QQ factor that sampleW_VY() has just built --
  // no refactorization, and no separate chain per group.
  //
  // Both quantities are free of the group's OWN measurement mixing variable:
  // the downdate removes exactly the term that carried it, so
  // QQ_-I = QQ - A_I' S_I^-1 A_I and b_-I = b - A_I' S_I^-1 r_I no longer
  // depend on noise_V_I. The caller therefore marginalises noise_V_I over its
  // prior when forming the predictive and the importance weight.
  //
  // `groups` holds 0-based observation indices. `chunk_cols` bounds the width
  // of the dense right-hand side handed to the solver at once. Results go into
  // plain C++ buffers -- out_mean[g] is k x n and out_cov[g] is the packed
  // lower triangle, k(k+1)/2 x n, both column-major by draw -- so that this can
  // run off the main thread, where no R object may be allocated.
  // The (W, V, noise_V) the model was constructed with -- for a fitted object
  // that is the state estimation left behind, which is far closer to each
  // fold's posterior than the prior is. Saved once and restored before every
  // fold, so folds start identically regardless of the order (or the thread)
  // they run in.
  void snapshot_state();
  void restore_state();

  // Drop observations (0-based) from the likelihood by zeroing their
  // measurement precision; clear_obs_mask() restores them.
  void set_obs_mask(const std::vector<int> &drop);
  void clear_obs_mask();
  // Seeds the block AND every latent's own stream. Seeding only the block
  // left each latent carrying state from whatever chain the worker ran before,
  // so a fold's draws depended on how the parallel loop was scheduled.
  void reseed(unsigned long seed) {
    rng.seed(seed);
    for (size_t i = 0; i < latents.size(); ++i)
      latents[i]->reseed(seed + 7919UL * (unsigned long)(i + 1));
  }

  // One exact leave-group-out chain: mask the group out, run the sampler, and
  // record eta_I = (AZ W)_I for every retained draw. This is what
  // cross_validation() does per fold, without leaving C++ or rebuilding the
  // model. `out_eta` is |I| x n, column-major by draw.
  // `start_W` overrides the snapshot's W for this chain. Passing a DIFFERENT
  // one per chain (the per-chain states the fit already stored) is what makes
  // the between-chain spread a real convergence signal: chains started from the
  // same point can agree while all being stuck.
  void loo_chain(const std::vector<int> &drop, int n, int n_burnin,
                 std::vector<double> &out_eta,
                 const VectorXd *start_W = nullptr);

  void group_cv_raw(const std::vector<std::vector<int>> &groups, int n,
                    int n_burnin, int chunk_cols, bool allow_inner_threads,
                    std::vector<std::vector<double>> &out_mean,
                    std::vector<std::vector<double>> &out_cov);

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
