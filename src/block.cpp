// Implementation for block model and block_rep

#include "block.h"
#include "include/phase_timing.h"
#include "include/factor_counters.h"
#include "include/solver.h"
#include "include/thread_io.h"
#include "prior.h"
#include "sample_rGIG.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>

using std::pow;

namespace {
// theta_sigma keeps every column of B_sigma, but the optimizer carries only the
// unfixed components and every gradient and Hessian entry is indexed by that
// shorter set. Anything assembled at the full width of B_sigma therefore has to
// be mapped through this list; taking the leading columns instead is correct
// only when the free components happen to come first.
std::vector<int> free_sigma_cols(const std::vector<bool> &fixed, int n_cols) {
  std::vector<int> free_cols;
  if ((int)fixed.size() == n_cols) {
    for (int i = 0; i < n_cols; ++i)
      if (!fixed[i])
        free_cols.push_back(i);
  } else {
    for (int i = 0; i < n_cols; ++i)
      free_cols.push_back(i);
  }
  return free_cols;
}
} // namespace

namespace ngme_counters {
std::atomic<long long> QQ_builds{0};
std::atomic<long long> probe_solves{0};
std::atomic<long long> gibbs_passes{0};
std::atomic<long long> fisher_solves{0};
std::atomic<long long> k_probe_solves{0};
thread_local probe_role current_probe_role = probe_role::qq;
std::atomic<long long> QQ_analyzes{0};
std::atomic<long long> K_analyzes{0};
void reset_all() {
  QQ_builds.store(0, std::memory_order_relaxed);
  probe_solves.store(0, std::memory_order_relaxed);
  gibbs_passes.store(0, std::memory_order_relaxed);
  fisher_solves.store(0, std::memory_order_relaxed);
  k_probe_solves.store(0, std::memory_order_relaxed);
  QQ_analyzes.store(0, std::memory_order_relaxed);
  K_analyzes.store(0, std::memory_order_relaxed);
}
namespace {
bool env_flag(const char *name) {
  const char *v = std::getenv(name);
  if (v == nullptr || v[0] == '\0')
    return false;
  const std::string s(v);
  return !(s == "0" || s == "false" || s == "FALSE");
}
} // namespace
bool cache_disabled() {
  // Re-read every call so that a test can toggle it with Sys.setenv() inside a
  // running session. getenv() costs nothing next to a sparse factorization.
  return env_flag("NGME2_DISABLE_FACTOR_CACHE");
}
} // namespace ngme_counters

namespace {
std::string parse_prior_target(const Rcpp::List &prior_list,
                               const std::string &default_target = "coef") {
  std::string target = prior_list.containsElementNamed("target")
                           ? Rcpp::as<std::string>(prior_list["target"])
                           : default_target;
  if (target != "coef" && target != "field") {
    throw std::invalid_argument("prior target must be 'coef' or 'field'");
  }
  return target;
}

VectorXd prior_score_vec(const std::string &type, const VectorXd &param,
                         const VectorXd &x) {
  VectorXd out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out(i) = PriorUtil::d_log_dens(type, param, x(i));
  }
  return out;
}

void parse_prior_spec(const Rcpp::List &prior_list, std::string &type,
                      VectorXd &param, std::string &target) {
  type = Rcpp::as<std::string>(prior_list["type"]);
  param = Rcpp::as<VectorXd>(prior_list["param"]);
  target = parse_prior_target(prior_list);
}
} // namespace

// -------------- Block Model class ----------------
BlockModel::BlockModel(const Rcpp::List &block_model, unsigned long seed)
    : rng(seed), X(Rcpp::as<MatrixXd>(block_model["X"])),
      Y(Rcpp::as<VectorXd>(block_model["Y"])),
      W_sizes(Rcpp::as<int>(block_model["W_sizes"])),
      V_sizes(Rcpp::as<int>(block_model["V_sizes"])),
      beta(Rcpp::as<VectorXd>(block_model["feff"])), n_obs(Y.size()),
      n_la_params(Rcpp::as<int>(block_model["n_la_params"])),
      n_feff(beta.size()), n_merr(Rcpp::as<int>(block_model["n_merr"])),
      n_repl(Rcpp::as<int>(block_model["n_repl"])), Q_eps(n_obs, n_obs),
      dQ_eps(n_obs, n_obs), n_params(n_la_params + n_feff + n_merr),

      debug(false), A(n_obs, W_sizes), K(V_sizes, W_sizes), Q(W_sizes, W_sizes),
      QQ(W_sizes, W_sizes),

      p_vec(n_obs), a_vec(n_obs), b_vec(n_obs), noise_V(VectorXd::Ones(n_obs)),
      noise_prevV(VectorXd::Ones(n_obs)),

      curr_iter(0), all_gaussian(Rcpp::as<bool>(block_model["all_gaussian"])),
      rao_blackwell(Rcpp::as<bool>(Rcpp::as<Rcpp::List>(
          block_model["control_ngme"])["rao_blackwellization"])),
      par_names(Rcpp::as<std::vector<std::string>>(block_model["par_names"])) {
  // 1. Init controls
  Rcpp::List control_ngme = block_model["control_ngme"];
  // const double stepsize = control_ngme["stepsize"];
  // bool init_sample_W = Rcpp::as<bool>(control_ngme["init_sample_W"]);
  n_gibbs = Rcpp::as<int>(control_ngme["n_gibbs_samples"]);
  int n_trace_iter = Rcpp::as<int>(control_ngme["n_trace_iter"]);
  // The operator traces are solved against K, not against the full block
  // precision, so they cost a fraction of a QQ probe while feeding both the
  // gradient of the operator parameters and the H_K block of the
  // preconditioner. Carrying them at a separate budget is what lets the spend
  // be split between the two rather than tied one-to-one.
  int n_trace_iter_k = n_trace_iter;
  if (control_ngme.containsElementNamed("n_trace_iter_k")) {
    Rcpp::IntegerVector v = control_ngme["n_trace_iter_k"];
    if (v.size() > 0 && v[0] != NA_INTEGER && v[0] > 0)
      n_trace_iter_k = v[0];
  }
  auto map_solver_type = [](int backend, int factor) {
    switch (backend) {
    case 0: // eigen
      return factor == 0 ? 0 : 1;
    case 1:                  // cholmod
      return factor == 0 ? 2 /*LLT (supernodal)*/
                         : 3 /*LDLT via CholmodDecomposition*/;
    case 2: // accelerate
      return factor == 0 ? 4 : 5;
    case 3: // pardiso
      return factor == 0 ? 6 : 7;
    default:
      throw std::invalid_argument("solver_backend out of range (expected 0-3)");
    }
  };

  int solver_backend = control_ngme.containsElementNamed("solver_backend")
                           ? Rcpp::as<int>(control_ngme["solver_backend"])
                           : 0;
  int solver_factor = control_ngme.containsElementNamed("solver_factor")
                          ? Rcpp::as<int>(control_ngme["solver_factor"])
                          : 0;
  int solver_type = map_solver_type(solver_backend, solver_factor);
  // 0 = LU of K, 1 = Cholesky of the normal equations, for non-symmetric
  // operators; see control_opt(nonsym_solver=).
  int nonsym_solver = control_ngme.containsElementNamed("nonsym_solver")
                          ? Rcpp::as<int>(control_ngme["nonsym_solver"])
                          : 0;
  robust = control_ngme.containsElementNamed("robust")
               ? Rcpp::as<bool>(control_ngme["robust"])
               : false;

  nig_param_std = control_ngme.containsElementNamed("nig_param_std")
                      ? Rcpp::as<int>(control_ngme["nig_param_std"])
                      : 0;
  precond_meas_sigma_ =
      control_ngme.containsElementNamed("precond_meas_sigma")
          ? Rcpp::as<int>(control_ngme["precond_meas_sigma"])
          : 0;
  fisher_refresh_every_ =
      control_ngme.containsElementNamed("fisher_refresh_every")
          ? std::max(1, Rcpp::as<int>(control_ngme["fisher_refresh_every"]))
          : 10;
  // reduce_var    =  Rcpp::as<bool>   (control_ngme["reduce_var"]);
  // reduce_power  =  Rcpp::as<double> (control_ngme["reduce_power"]);
  // threshold   =  Rcpp::as<double> (control_ngme["threshold"]);

  if (debug)
    ngme_io::out() << "Begin Block Constructor" << std::endl;

  // 2. Init Fixed effects
  bool fix_beta = control_ngme.containsElementNamed("fix_beta")
                      ? Rcpp::as<bool>(control_ngme["fix_beta"])
                      : Rcpp::as<bool>(control_ngme["fix_feff"]);
  fix_flag[block_fix_beta] = fix_beta;
  if (beta.size() == 0)
    fix_flag[block_fix_beta] = true;

  // init priors for fixed effects
  prior_beta_type.resize(n_feff, "none");
  prior_beta_param.resize(n_feff, VectorXd::Zero(0));
  prior_beta_target.resize(n_feff, "coef");
  if (block_model.containsElementNamed("prior_beta")) {
    Rcpp::List prior_beta_list = Rcpp::as<Rcpp::List>(block_model["prior_beta"]);
    if (n_feff > 0 && prior_beta_list.size() > 0) {
      if (prior_beta_list.containsElementNamed("type")) {
        std::string t;
        VectorXd p;
        std::string target;
        parse_prior_spec(prior_beta_list, t, p, target);
        for (int i = 0; i < n_feff; ++i) {
          prior_beta_type[i] = t;
          prior_beta_param[i] = p;
          prior_beta_target[i] = target;
        }
      } else {
        if (prior_beta_list.size() != n_feff) {
          throw std::invalid_argument(
              "prior_beta must have length equal to number of fixed effects");
        }
        for (int i = 0; i < n_feff; ++i) {
          Rcpp::List one = Rcpp::as<Rcpp::List>(prior_beta_list[i]);
          parse_prior_spec(one, prior_beta_type[i], prior_beta_param[i],
                           prior_beta_target[i]);
        }
      }
    }
  }
  for (int i = 0; i < n_feff; ++i) {
    if (prior_beta_target[i] != "coef") {
      throw std::invalid_argument(
          "beta prior target='field' is not supported; use target='coef'");
    }
  }

  // 4. Init latent models
  Rcpp::List latents_in = block_model["models"];
  n_latent = latents_in.size(); // how many latent model
  if (n_latent == 0)
    rao_blackwell = false;
  for (int i = 0; i < n_latent; ++i) {
    // construct acoording to models
    Rcpp::List latent_in = Rcpp::as<Rcpp::List>(latents_in[i]);
    latent_in["solver_type"] = solver_type;
    latent_in["nonsym_solver"] = nonsym_solver;
    latent_in["n_trace_iter"] = n_trace_iter_k;
    // The operator traces use the same fill gate as the QQ traces do, so the
    // threshold has to reach Operator::update_all as well. Read from
    // control_ngme here rather than from the member, which is assigned further
    // down in this constructor.
    // The operator traces are taken once per optimizer iteration, while the
    // block-precision traces are retaken on every Gibbs pass. The exact path
    // therefore costs the operator far less over a run, so the two gates are
    // separate rather than sharing one threshold.
    latent_in["selinv_max_fill"] =
        control_ngme.containsElementNamed("selinv_max_fill_k")
            ? Rcpp::as<double>(control_ngme["selinv_max_fill_k"])
        : control_ngme.containsElementNamed("selinv_max_fill")
            ? Rcpp::as<double>(control_ngme["selinv_max_fill"])
            : 4.0;
    latent_in["robust"] = robust;
    latent_in["selinv_cost_ratio"] =
        control_ngme.containsElementNamed("selinv_cost_ratio")
            ? Rcpp::as<double>(control_ngme["selinv_cost_ratio"])
            : 2.0;
    latent_in["trace_probing"] =
        control_ngme.containsElementNamed("trace_probing")
            ? Rcpp::as<bool>(control_ngme["trace_probing"])
            : true;
    latent_in["trace_probing_max_dist"] =
        control_ngme.containsElementNamed("trace_probing_max_dist")
            ? Rcpp::as<int>(control_ngme["trace_probing_max_dist"])
            : 4;
    // The operator budget follows the QQ budget when trace_adapt_k is on, so
    // the colouring search is capped at the same ceiling the QQ side uses.
    latent_in["trace_probing_max_colours"] =
        control_ngme.containsElementNamed("trace_adapt_max")
            ? Rcpp::as<int>(control_ngme["trace_adapt_max"])
            : n_trace_iter_k;
    latent_in["trace_probing_raise_budget"] =
        control_ngme.containsElementNamed("trace_probing_raise_budget")
            ? Rcpp::as<double>(control_ngme["trace_probing_raise_budget"])
            : 1.0;
    latent_in["nig_param_std"] = nig_param_std;
    unsigned long latent_seed = seed + (i + 1) * 1000;
    latents.push_back(std::make_shared<Latent>(latent_in, latent_seed));
  }

  if (debug)
    ngme_io::out() << "before set block A" << std::endl;
  /* Init A */
  int n = 0;
  for (std::vector<std::shared_ptr<Latent>>::iterator it = latents.begin();
       it != latents.end(); it++) {
    setSparseBlock(&A, 0, n, (*it)->getA());
    n += (*it)->get_W_size();
  }
  if (debug)
    ngme_io::out() << "After set block K" << std::endl;

  // 5. Init measurement noise (consider corr_measure)
  Rcpp::List noise_in = block_model["noise"];
  fix_flag[block_fix_theta_mu] = Rcpp::as<bool>(noise_in["fix_theta_mu"]);

  // Handle vector-based fix_theta_sigma - simplified for block model
  Rcpp::LogicalVector fix_theta_sigma_r =
      Rcpp::as<Rcpp::LogicalVector>(noise_in["fix_theta_sigma"]);
  fix_theta_sigma_vec =
      std::vector<bool>(fix_theta_sigma_r.begin(), fix_theta_sigma_r.end());
  fix_flag[block_fix_theta_sigma] =
      std::all_of(fix_theta_sigma_r.begin(), fix_theta_sigma_r.end(),
                  [](bool x) { return x; });

  fix_flag[blcok_fix_V] = Rcpp::as<bool>(noise_in["fix_V"]);
  fix_flag[block_fix_theta_nu] = noise_in.containsElementNamed("fix_theta_nu")
                                     ? Rcpp::as<bool>(noise_in["fix_theta_nu"])
                                     : false;
  fix_flag[block_fix_rho] = noise_in.containsElementNamed("fix_rho")
                                ? Rcpp::as<bool>(noise_in["fix_rho"])
                                : false;
  // shared_sigma has been removed; always use standard scaling with sigma^2 and
  // V

  B_mu = (Rcpp::as<MatrixXd>(noise_in["B_mu"])),
  theta_mu = (Rcpp::as<VectorXd>(noise_in["theta_mu"])),
  n_theta_mu = (Rcpp::as<int>(noise_in["n_theta_mu"])),

  B_sigma = (Rcpp::as<MatrixXd>(noise_in["B_sigma"])),
  theta_sigma = (Rcpp::as<VectorXd>(noise_in["theta_sigma"])),
  n_theta_sigma = (Rcpp::as<int>(noise_in["n_theta_sigma"])),

  B_nu = (Rcpp::as<MatrixXd>(noise_in["B_nu"])),
  theta_nu = (Rcpp::as<VectorXd>(noise_in["theta_nu"])),
  n_theta_nu = (Rcpp::as<int>(noise_in["n_theta_nu"])),

  rb_trace_noise_sigma = VectorXd::Zero(n_theta_sigma),

  nu_lower_bound = noise_in.containsElementNamed("nu_lower_bound")
                       ? Rcpp::as<double>(noise_in["nu_lower_bound"])
                       : 0.0,

  family = Rcpp::as<string>(noise_in["noise_type"]);
  {
    Rcpp::List cn = block_model["control_ngme"];
    if (cn.containsElementNamed("debug"))
      debug = Rcpp::as<bool>(cn["debug"]);
#ifdef __APPLE__
    // Fill-reducing ordering for the Accelerate backend; see solver.h. Process-
    // global because it only has to reach the solvers this fit constructs, and
    // it does not change during a fit.
    if (cn.containsElementNamed("solver_order"))
      ngme_set_accel_order(Rcpp::as<int>(cn["solver_order"]));
#endif
    if (cn.containsElementNamed("selinv_max_fill"))
      selinv_max_fill = Rcpp::as<double>(cn["selinv_max_fill"]);
    if (cn.containsElementNamed("n_fisher_probes")) {
      Rcpp::IntegerVector v = cn["n_fisher_probes"];
      if (v.size() > 0 && v[0] != NA_INTEGER && v[0] > 0)
        n_fisher_probes_ = v[0];
    }
    if (cn.containsElementNamed("trace_adapt"))
      trace_adapt = Rcpp::as<bool>(cn["trace_adapt"]);
    if (cn.containsElementNamed("trace_adapt_frac"))
      trace_adapt_frac = Rcpp::as<double>(cn["trace_adapt_frac"]);
    if (cn.containsElementNamed("trace_adapt_min"))
      trace_adapt_min = Rcpp::as<int>(cn["trace_adapt_min"]);
    if (cn.containsElementNamed("trace_adapt_max"))
      trace_adapt_max = Rcpp::as<int>(cn["trace_adapt_max"]);
    if (cn.containsElementNamed("trace_probing"))
      trace_probing = Rcpp::as<bool>(cn["trace_probing"]);
    if (cn.containsElementNamed("trace_probing_max_dist"))
      trace_probing_max_dist = Rcpp::as<int>(cn["trace_probing_max_dist"]);
    if (cn.containsElementNamed("trace_adapt_rule"))
      trace_adapt_cost_rule =
          Rcpp::as<std::string>(cn["trace_adapt_rule"]) == "cost";
    if (cn.containsElementNamed("trace_probing_raise_budget"))
      trace_probing_raise_budget =
          Rcpp::as<double>(cn["trace_probing_raise_budget"]);
    if (cn.containsElementNamed("trace_adapt_k"))
      trace_adapt_k = Rcpp::as<bool>(cn["trace_adapt_k"]);
  }

  noise_mu = B_mu * theta_mu;
  noise_sigma = (B_sigma * theta_sigma).array().exp();
  noise_nu = nu_lower_bound + (B_nu * theta_nu).array().exp();

  rho = Rcpp::as<VectorXd>(noise_in["rho"]);
  n_rho = noise_in.containsElementNamed("n_rho")
              ? Rcpp::as<int>(noise_in["n_rho"])
              : 0;
  corr_measure = Rcpp::as<bool>(noise_in["corr_measurement"]);

  // init priors for noise_parameter
  Rcpp::List prior_list = Rcpp::as<Rcpp::List>(noise_in["prior_mu"]);
  parse_prior_spec(prior_list, prior_mu_type, prior_mu_param, prior_mu_target);
  prior_list = Rcpp::as<Rcpp::List>(noise_in["prior_sigma"]);
  parse_prior_spec(prior_list, prior_sigma_type, prior_sigma_param,
                   prior_sigma_target);
  prior_list = Rcpp::as<Rcpp::List>(noise_in["prior_nu"]);
  parse_prior_spec(prior_list, prior_nu_type, prior_nu_param, prior_nu_target);

  if (family != "normal") {
    NoiseUtil::update_gig(family, noise_nu, p_vec, a_vec, b_vec);
  }

  if (corr_measure) {
    cor_cols = Rcpp::as<vector<int>>(noise_in["cor_cols"]);
    cor_rows = Rcpp::as<vector<int>>(noise_in["cor_rows"]);
    has_correlation = Rcpp::as<vector<bool>>(noise_in["has_correlation"]);

    n_corr_pairs = Rcpp::as<int>(noise_in["n_corr_pairs"]);
    vector<Triplet<double>> Q_eps_triplet, dQ_eps_triplet;
    for (int i = 0; i < cor_cols.size(); ++i) {
      Q_eps_triplet.push_back(Triplet<double>(cor_rows[i], cor_cols[i],
                                              cor_rows[i] == cor_cols[i]));
      if (has_correlation[cor_rows[i]]) {
        // ignore uncorrelated locations
        dQ_eps_triplet.push_back(Triplet<double>(cor_rows[i], cor_cols[i],
                                                 cor_rows[i] == cor_cols[i]));
      }
    }
    SparseMatrix<double> Q_eps_lower(n_obs, n_obs);
    Q_eps_lower.setFromTriplets(Q_eps_triplet.begin(), Q_eps_triplet.end());
    Q_eps = Q_eps_lower.selfadjointView<Lower>();
    sqrt_Rinv = Q_eps; // initialize sqrt_Rinv as Q_eps

    SparseMatrix<double> dQ_lower(n_obs, n_obs);
    dQ_lower.setFromTriplets(dQ_eps_triplet.begin(), dQ_eps_triplet.end());
    dQ_eps = dQ_lower.selfadjointView<Lower>();

    update_Q_eps(rho(0));
    // std::cout << "rho(0) = " << rho(0) << std::endl;
    // std::cout << "noise_sigma = " << noise_sigma.transpose() << std::endl;
    // std::cout << "Init Q_eps: \n" << Q_eps << std::endl;
  }

  if (debug)
    ngme_io::out() << "After block construct noise" << std::endl;

  // 6. Fix V and init V
  if (noise_in.containsElementNamed("V") && !Rf_isNull(noise_in["V"])) {
    noise_V = Rcpp::as<VectorXd>(noise_in["V"]);
    noise_prevV = noise_V;
  }

  // 7. Init solvers
  assemble();
  if (debug)
    ngme_io::out() << "After assemble" << std::endl;

  if (n_latent > 0) {
    VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
    Q = K.transpose() * inv_SV.asDiagonal() * K;
    // Build AZ for measurement term
    const SparseMatrix<double> &AZ = get_AZ();
    if (!corr_measure) {
      QQ = Q + AZ.transpose() *
                   noise_sigma.array()
                       .pow(-2)
                       .matrix()
                       .cwiseQuotient(noise_V)
                       .asDiagonal() *
                   AZ;
    } else {
      QQ = Q + AZ.transpose() * Q_eps * AZ;
    }

    // Initialize solver with requested backend and Hutchinson iters; QQ is SPD
    chol_QQ.init(QQ.rows(), n_trace_iter, /*symmetric*/ true, solver_type);
    chol_QQ.analyze(QQ);
    setup_qq_probing();
    record_QQ_pattern();
    chol_QQ.compute(QQ);
    // Deliberately leave QQ_valid == false: this initial assembly is not
    // identical to update_QQ(), which additionally clamps inv_SV from below
    // and, under `robust`, symmetrizes and jitters QQ. The first sampleW_VY()
    // therefore still does a full update_QQ(), exactly as before caching; the
    // recorded pattern only spares that call a redundant symbolic phase.
  }

  // 8. optimizer related
  // stepsizes = VectorXd::Constant(n_params, stepsize);
  steps_to_threshold = VectorXd::Constant(n_params, 0);
  indicate_threshold = VectorXd::Constant(n_params, 0);

  // if (debug) std::cout << "After init solver && before sampleW_V" <<
  // std::endl;

  // Default: enable RB for all-Gaussian models and use conditional W
  if (rao_blackwell) {
    // Initialize block_dKs of length n_latent
    block_dK.resize(n_latent);
    for (int i = 0; i < n_latent; i++) {
      block_dK[i].resize(latents[i]->get_n_theta_K());
      for (int j = 0; j < latents[i]->get_n_theta_K(); j++) {
        block_dK[i][j] = SparseMatrix<double>(V_sizes, W_sizes);
      }
    }
  }

  // Initialize gradient covariance matrix
  grad_covariance = MatrixXd::Zero(n_params, n_params);

  if (debug)
    ngme_io::out() << "End Block Constructor" << std::endl;
}

void BlockModel::burn_in(int iterations) {
  if (debug)
    ngme_io::out() << "Start burn-in for " << iterations
                << " iterations of burn-in" << std::endl;
  for (int i = 0; i < iterations; i++) {
    sample_cond_V();
    sampleW_VY(true);
    sample_cond_noise_V();
  }
  if (debug)
    ngme_io::out() << "End burn-in" << std::endl;
}

void BlockModel::setW(const VectorXd &W) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it =
           latents.begin();
       it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->setW(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::set_cond_W(const VectorXd &W) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it =
           latents.begin();
       it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->set_cond_W(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::setPrevW(const VectorXd &W) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it =
           latents.begin();
       it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->setPrevW(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::setPrevV(const VectorXd &V) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it =
           latents.begin();
       it != latents.end(); it++) {
    int size = (*it)->get_V_size();
    (*it)->setPrevV(V.segment(pos, size));
    pos += size;
  }
}

// sample W|VY
void BlockModel::snapshot_state() {
  snap_W.clear(); snap_V.clear();
  snap_W.reserve(n_latent); snap_V.reserve(n_latent);
  for (int li = 0; li < n_latent; ++li) {
    snap_W.push_back(latents[li]->getW());
    snap_V.push_back(latents[li]->getV());
  }
  snap_noise_V = noise_V;
  snap_taken = true;
}

void BlockModel::restore_state() {
  if (!snap_taken)
    return;
  for (int li = 0; li < n_latent; ++li) {
    latents[li]->setW(snap_W[li]);
    latents[li]->setV(snap_V[li]);
  }
  noise_V = snap_noise_V;
  invalidate_QQ();
  invalidate_measurement();
}

// ---- exact leave-group-out support ----
//
// Dropping observations by zeroing their measurement precision is exact: the
// likelihood contribution of observation j enters QQ only through H'H with
// H = sqrt(D) AZ, and M only through (AZ)' D r, so D_j = 0 removes it from
// both. The latent block, the operator and the symbolic factorisation are
// untouched, and removing rows can only shrink QQ's pattern, never grow it.
void BlockModel::set_obs_mask(const std::vector<int> &drop) {
  obs_weight = VectorXd::Ones(n_obs);
  for (int idx : drop) {
    if (idx < 0 || idx >= n_obs)
      throw std::runtime_error("set_obs_mask(): index out of range");
    obs_weight(idx) = 0.0;
  }
  invalidate_measurement();
}

void BlockModel::clear_obs_mask() {
  if (obs_weight.size() == n_obs && obs_weight.isApprox(VectorXd::Ones(n_obs)))
    return;
  obs_weight = VectorXd::Ones(n_obs);
  invalidate_measurement();
}

// One exact leave-group-out chain, equivalent to what cross_validation() runs
// per fold but without leaving C++ or rebuilding the model.
void BlockModel::loo_chain(const std::vector<int> &drop, int n, int n_burnin,
                           std::vector<double> &out_eta,
                           const VectorXd *start_W) {
  const int k = (int)drop.size();
  out_eta.assign((size_t)k * n, 0.0);
  if (n_latent == 0)
    return;

  set_obs_mask(drop);
  // Start every fold from the state the model was BUILT with -- for a fitted
  // object, the W and V estimation left behind. Resetting to the prior instead
  // (W = 0, V drawn from its prior) was measurably worse: it is far from the
  // fold's posterior, so it needs much more burn-in, and it is a worse start
  // than cross_validation() gets, since a freshly constructed BlockModel keeps
  // the fitted W. Restoring a snapshot keeps folds identical and independent of
  // the order they run in, which inheriting the previous fold's state would not.
  restore_state();
  if (start_W != nullptr && start_W->size() == W_sizes) {
    setW(*start_W);
    invalidate_QQ();
  }

  rao_blackwell = false;
  if (!all_gaussian)
    burn_in(n_burnin);
  else
    for (int i = 0; i < n_burnin; ++i)
      sampleW_VY(true);

  // AZ is column-major, so AZ.row(i) has to walk every column: using it inside
  // the draw loop costs O(nnz(AZ)) per held-out observation per draw instead of
  // O(nnz of that row). Extract the k rows we need with ONE pass over AZ, then
  // the per-draw work is k short dot products. (group_cv_raw() builds the same
  // index for the same reason.)
  const SparseMatrix<double> &AZ = get_AZ();
  std::vector<int> row_of(n_obs, -1);
  for (int i = 0; i < k; ++i)
    row_of[drop[i]] = i;
  std::vector<std::vector<std::pair<int, double>>> rows(k);
  for (int c = 0; c < AZ.outerSize(); ++c)
    for (SparseMatrix<double>::InnerIterator it(AZ, c); it; ++it) {
      const int r = row_of[it.row()];
      if (r >= 0)
        rows[r].emplace_back((int)c, it.value());
    }

  for (int it = 0; it < n; ++it) {
    if (!all_gaussian)
      sample_cond_V();
    sampleW_VY();
    sample_cond_noise_V(true);
    const VectorXd W = getW();
    const double *wp = W.data();
    for (int i = 0; i < k; ++i) {
      double acc = 0.0;
      for (const auto &e : rows[i])
        acc += e.second * wp[e.first];
      out_eta[(size_t)k * it + i] = acc;
    }
  }
  clear_obs_mask();
}

void BlockModel::sampleW_VY(bool burn_in) {
  ngme_timing::Scope _sw(ngme_timing::samplew_us());
  if (n_latent == 0)
    return;
  // Ensure QQ is consistent with current K, Z, and measurement precision before
  // sampling. Nothing that enters QQ changes between the Gibbs draws of a
  // Gaussian model, so this is a no-op after the first draw of an iteration.
  { ngme_timing::Scope _s(ngme_timing::sw_ensureQQ_us()); ensure_QQ(); }
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());

  // M = K' * inv(SV) * mean + Z'^ A'^ inv(Sigma) * (Y - X * beta - (1 - V) mu)
  ngme_timing::Scope _sm(ngme_timing::sw_M_us());
  // Parenthesised so the DIAGONAL meets the vector first.
  //
  // `A.transpose() * d.asDiagonal() * v` groups as `(A^T * d) * v`, so Eigen
  // builds a scaled transpose of A. Scaling the vector instead leaves a single
  // sparse-transpose-times-vector, which Eigen does without materialising anything.
  VectorXd M = K.transpose() * (inv_SV.asDiagonal() * getMean());

  const SparseMatrix<double> &AZ = get_AZ();
  if (!corr_measure) {
    M += AZ.transpose() *
         (meas_prec().asDiagonal() * get_residual_part());
  } else {
    M += AZ.transpose() * (Q_eps * get_residual_part());
  }

  _sm.stop();
  const SparseMatrix<double> *Gp, *Hp;
  { ngme_timing::Scope _s(ngme_timing::sw_G_us()); Gp = &get_G(inv_SV); }
  { ngme_timing::Scope _s(ngme_timing::sw_H_us()); Hp = &get_sqrt_AtSVA(); }
  const SparseMatrix<double> &G = *Gp;
  const SparseMatrix<double> &H = *Hp;
  unsigned long seed1 = rng();
  unsigned long seed2 = rng();
  VectorXd z1 = NoiseUtil::rnorm_vec(G.rows(), 0, 1, seed1);
  VectorXd z2 = NoiseUtil::rnorm_vec(H.rows(), 0, 1, seed2);

  // sample W ~ N(QQ^-1*M, QQ^-1)

  // Sampling method using cholesky decomposition, using matrixL()
  // W = QQ^-1*M + L^-1*z1 + L^-T*z2
  // VectorXd W = chol_QQ.rMVN(M, z1);

  // Sampling method using tricks, purely solve, not matrixL()
  { ngme_timing::Scope _s(ngme_timing::rmvn_us());
    VectorXd W_draw = chol_QQ.rMVN(G, H, M, z1, z2);
    setW(W_draw);

    if (rao_blackwell && !burn_in) {
      VectorXd W_mean = chol_QQ.solve(M);
      set_cond_W(W_mean);
    } }
  // std::cout << "size of W and time of sampling is " << W.size() << " " <<
  // time << std::endl; if (debug) std::cout << "Finish sampling W" <<
  // std::endl;
}

// ---------------- get, set update gradient ------------------
// order is Latent, merr, feff
VectorXd BlockModel::get_parameter() {
  if (debug)
    ngme_io::out() << "Start block get_parameter" << std::endl;
  VectorXd thetas(n_params);
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it =
           latents.begin();
       it != latents.end(); it++) {
    VectorXd theta = (*it)->get_parameter();
    thetas.segment(pos, theta.size()) = theta;
    pos += theta.size();
  }
  thetas.segment(n_la_params, n_merr) = get_theta_merr();

  if (!fix_flag[block_fix_beta])
    thetas.segment(n_la_params + n_merr, n_feff) = beta;

  if (debug)
    ngme_io::out() << "Finish block get_parameter" << std::endl;
  return thetas;
}

void BlockModel::set_parameter_and_update(const VectorXd &Theta,
                                          bool with_precond) {
  if (debug)
    ngme_io::out() << "Start block set_parameter" << std::endl;
  // std::chrono::steady_clock::time_point startTime, endTime; startTime =
  // std::chrono::steady_clock::now();
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::iterator it = latents.begin();
       it != latents.end(); it++) {
    int theta_len = (*it)->get_n_params();
    VectorXd theta = Theta.segment(pos, theta_len);
    (*it)->set_parameter_and_update(theta, with_precond);
    pos += theta_len;
  }

  // measurement noise
  set_theta_merr(Theta.segment(n_la_params, n_merr));
  if (family != "normal")
    NoiseUtil::update_gig(family, noise_nu, p_vec, a_vec, b_vec);

  // fixed effects
  if (!fix_flag[block_fix_beta]) {
    beta = Theta.segment(n_la_params + n_merr, n_feff);
  }

  assemble(); // update K,dK,d2K after
  // K, Z, the latent sigmas and the measurement precision may all have moved,
  // so both the AZ block matrix and QQ have to be rebuilt on next use.
  invalidate_AZ();
  // endTime = std::chrono::steady_clock::now(); update_time =
  // std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
  // std::cout << "block set time (ms): " <<
  // std::chrono::duration_cast<std::chrono::milliseconds>(endTime -
  // startTime).count() << std::endl;
  curr_iter++;
  // QQ is refreshed lazily when sampling W in sampleW_VY().
}

// --------- Fiexed effects and Measurement Error ---------------
VectorXd BlockModel::grad_beta() {
  // Measurement precision: 1/(sigma^2 V), or Q_eps when correlated -- the same
  // weighting the W sampler and the field score s_full use.
  VectorXd residual = get_residual(rao_blackwell); // + X * beta;
  VectorXd grads;
  if (!corr_measure) {
    VectorXd noise_inv_SV =
        noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
    grads = X.transpose() * (noise_inv_SV.asDiagonal() * residual);
  } else {
    grads = X.transpose() * (Q_eps * residual);
  }

  for (int l = 0; l < n_feff; ++l) {
    grads(l) +=
        PriorUtil::d_log_dens(prior_beta_type[l], prior_beta_param[l], beta(l));
  }

  // shared_sigma removed: keep generic form only
  //  * residual.cwiseQuotient(noise_sigma);

  return grads;
}

VectorXd BlockModel::grad_theta_mu() {
  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());
  VectorXd residual = get_residual(rao_blackwell);
  VectorXd grad = VectorXd::Zero(n_theta_mu);
  for (int l = 0; l < n_theta_mu; l++) {
    grad(l) = (noise_V - VectorXd::Ones(n_obs))
                  .cwiseProduct(B_mu.col(l).cwiseQuotient(noise_SV))
                  .dot(residual);
  }

  if (prior_mu_target == "coef") {
    for (int l = 0; l < n_theta_mu; l++) {
      grad(l) += PriorUtil::d_log_dens(prior_mu_type, prior_mu_param,
                                       theta_mu(l));
    }
  } else {
    grad += B_mu.transpose() *
            prior_score_vec(prior_mu_type, prior_mu_param, noise_mu);
  }
  return grad;
}

VectorXd BlockModel::grad_theta_sigma() {
  VectorXd grad = VectorXd::Zero(n_theta_sigma);
  VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);

  VectorXd residual = get_residual(rao_blackwell);
  VectorXd vsq = (residual).array().pow(2).matrix().cwiseQuotient(noise_SV);
  VectorXd tmp1 = vsq - VectorXd::Ones(n_obs);
  VectorXd full_grad = B_sigma.transpose() * tmp1; // dℓ/dθ for generic case
  int j = 0;
  for (int i = 0; i < full_grad.size(); ++i) {
    if (!fix_theta_sigma_vec[i]) {
      grad(j++) = full_grad(i);
    }
  }

  // Rao-Blackwell term from log|QQ|: keep sign consistent with theta_K path
  // theta_K: grad_accum -= rb_trace_K; then return -grad_accum
  // Here: do grad -= rb_trace_noise_sigma; then return -grad
  if (rao_blackwell)
    grad += rb_trace_noise_sigma;

  if (prior_sigma_target == "coef") {
    j = 0;
    for (int l = 0; l < theta_sigma.size(); l++) {
      if (!fix_theta_sigma_vec[l]) {
        grad(j++) += PriorUtil::d_log_dens(prior_sigma_type, prior_sigma_param,
                                           theta_sigma(l));
      }
    }
  } else {
    VectorXd sigma_lp = B_sigma * theta_sigma;
    VectorXd full_prior =
        B_sigma.transpose() * prior_score_vec(prior_sigma_type, prior_sigma_param, sigma_lp);
    j = 0;
    for (int l = 0; l < full_prior.size(); l++) {
      if (!fix_theta_sigma_vec[l]) {
        grad(j++) += full_prior(l);
      }
    }
  }

  // Return gradient of negative log-likelihood
  return grad;
}

// Expected (Fisher) information for the free measurement theta_sigma, of the
// marginal likelihood with W integrated out, given V. With S = diag(sigma^2 V)
// and dS/dtheta_l = 2 diag(b_l) S,
//   F_lk = 1/2 tr(Sigma_Y^-1 dS_l Sigma_Y^-1 dS_k) = 2 tr(N B_l N B_k),
//   N = S^1/2 Sigma_Y^-1 S^1/2 = I - H QQ^-1 H^T,   H = D^1/2 A Z.
// The complete-data Hessian is ~2n whatever sigma is, while the marginal
// information vanishes like sigma^4 as sigma -> 0 and the score like sigma^2.
// Preconditioning with the former freezes a chain that has drifted onto the
// sigma -> 0 plateau; with the latter the step grows there and it escapes.
// N is symmetric with eigenvalues in (0, 1], so the Hutchinson estimate
// (Nz)^T B_l (N B_k z) stays accurate relative to F even where F is tiny --
// unlike expanding it into n - 2 tr(DP) + tr(DPDP), a difference of ~n terms.
MatrixXd BlockModel::fisher_theta_sigma() {
  std::vector<int> cols;
  for (int i = 0; i < (int)theta_sigma.size(); ++i)
    if (!fix_theta_sigma_vec[i])
      cols.push_back(i);
  const int p = (int)cols.size();
  MatrixXd F = MatrixXd::Zero(p, p);
  if (p == 0)
    return F;
  const bool masked = obs_weight.size() == n_obs;

  if (n_latent == 0 || !QQ_valid) {
    // No field: N = I.
    for (int l = 0; l < p; ++l)
      for (int k = 0; k <= l; ++k) {
        VectorXd bb = B_sigma.col(cols[l]).cwiseProduct(B_sigma.col(cols[k]));
        if (masked)
          bb = bb.cwiseProduct(obs_weight);
        F(l, k) = F(k, l) = 2.0 * bb.sum();
      }
    return F;
  }

  const SparseMatrix<double> &H = get_sqrt_AtSVA();
  const int N = std::max(1, n_fisher_probes_ > 0 ? n_fisher_probes_
                                                : chol_QQ.get_N_iter());
  std::mt19937 gen(static_cast<unsigned>(rng()));
  std::bernoulli_distribution coin(0.5);
  // Dropped observations get a zero probe entry, which removes them exactly.
  MatrixXd Zp(n_obs, N);
  for (int c = 0; c < N; ++c)
    for (int r = 0; r < n_obs; ++r)
      Zp(r, c) = (masked && obs_weight(r) == 0.0) ? 0.0
                                                  : (coin(gen) ? 1.0 : -1.0);
  auto applyN = [&](const MatrixXd &Y) {
    MatrixXd rhs = H.transpose() * Y;
    // Counted separately from the trace probes: same unit of work, but a
    // different consumer, and lumping them together hides whether the traces
    // are still probing.
    ngme_counters::add(ngme_counters::fisher_solves, N);
    MatrixXd s = chol_QQ.solve(rhs);
    return MatrixXd(Y - H * s);
  };
  const MatrixXd NZ = applyN(Zp);
  std::vector<MatrixXd> NBZ(p);
  for (int k = 0; k < p; ++k) {
    const VectorXd bk = B_sigma.col(cols[k]);
    NBZ[k] = bk.isOnes() ? NZ : applyN(bk.asDiagonal() * Zp);
  }
  for (int l = 0; l < p; ++l) {
    const VectorXd bl = B_sigma.col(cols[l]);
    for (int k = 0; k <= l; ++k) {
      double acc = (NZ.array() * (bl.asDiagonal() * NBZ[k]).array()).sum();
      F(l, k) = F(k, l) = 2.0 * acc / N;
    }
  }
  return F;
}

VectorXd BlockModel::get_theta_merr() const {
  VectorXd theta_merr = VectorXd::Zero(n_merr);

  if (!fix_flag[block_fix_theta_mu])
    theta_merr.segment(0, n_theta_mu) = theta_mu;
  if (!fix_flag[block_fix_theta_sigma]) {
    int pos = 0;
    for (int i = 0; i < theta_sigma.size(); ++i) {
      if (!fix_theta_sigma_vec[i]) {
        theta_merr(n_theta_mu + pos) = theta_sigma(i);
        pos += 1;
      }
    }
  }
  if (!fix_flag[block_fix_theta_nu])
    theta_merr.segment(n_theta_mu + n_theta_sigma, n_theta_nu) = theta_nu;

  if (corr_measure && !fix_flag[block_fix_rho])
    theta_merr(n_merr - 1) = rho2th(rho(0));

  // The optimiser sees the standardised coordinates when they apply.
  if (int mode = merr_nig_mode()) {
    Eigen::Vector3d native(theta_mu(0), theta_sigma(0), theta_nu(0));
    theta_merr.head(3) = nig_std::from_native(mode, native);
  }
  return theta_merr;
}

// The standardised map mixes mu, sigma and nu, so it needs all three free,
// scalar and stationary, nu unshifted, and uncorrelated noise.
int BlockModel::merr_nig_mode() const {
  if (nig_param_std == 0 || family != "nig" || corr_measure)
    return 0;
  if (n_theta_mu != 1 || n_theta_sigma != 1 || n_theta_nu != 1 ||
      theta_sigma.size() != 1 || nu_lower_bound != 0.0)
    return 0;
  if (fix_flag[block_fix_theta_mu] || fix_flag[block_fix_theta_sigma] ||
      fix_flag[block_fix_theta_nu])
    return 0;
  if (!B_mu.isOnes() || !B_sigma.isOnes() || !B_nu.isOnes())
    return 0;
  return nig_param_std;
}

VectorXd BlockModel::grad_theta_merr() {
  VectorXd grad = VectorXd::Zero(n_merr);

  if (!fix_flag[block_fix_theta_mu])
    grad.segment(0, n_theta_mu) = grad_theta_mu();
  if (!fix_flag[block_fix_theta_sigma])
    grad.segment(n_theta_mu, n_theta_sigma) = grad_theta_sigma();
  if (!fix_flag[block_fix_theta_nu]) {
    grad.segment(n_theta_mu + n_theta_sigma, n_theta_nu) =
        -NoiseUtil::grad_theta_nu(family, B_nu, noise_nu, noise_V, noise_prevV,
                                  VectorXd::Ones(noise_V.size()),
                                  nu_lower_bound);
    VectorXd prior_grad = VectorXd::Zero(n_theta_nu);
    if (prior_nu_target == "coef") {
      for (int l = 0; l < n_theta_nu; ++l) {
        prior_grad(l) += PriorUtil::d_log_dens(prior_nu_type, prior_nu_param,
                                               theta_nu(l));
      }
    } else {
      VectorXd nu_lp = B_nu * theta_nu;
      prior_grad +=
          B_nu.transpose() * prior_score_vec(prior_nu_type, prior_nu_param, nu_lp);
    }
    grad.segment(n_theta_mu + n_theta_sigma, n_theta_nu) += prior_grad;
  }

  // grad of theta_rho
  if (corr_measure && !fix_flag[block_fix_rho]) {
    // Q_eps_solver.factorize(Q_eps);
    double trace = 0.5 * 2 * rho(0) / (1 - rho(0) * rho(0)) * n_corr_pairs;
    VectorXd res = get_residual();
    double drhs = -0.5 * (res).dot(dQ_eps * res);
    grad(n_merr - 1) = trace + drhs;
    grad(n_merr - 1) *= dtheta_th(rho(0));
    // std::cout << "drhs = " << drhs << std::endl;
    // std::cout << "trace = " << trace << std::endl;
    // std::cout << "grad of rho=" << grad(n_merr-1) << std::endl;
  }

  return grad;
}

void BlockModel::set_theta_merr(const VectorXd &theta_merr_in) {
  // if (debug) std::cout << "start set theta_merr" << std::endl;
  // noise_sigma / Q_eps feed the measurement block of QQ.
  invalidate_measurement();
  // The optimiser may be working in the standardised coordinates; everything
  // below (and every derivative) stays native, so convert on the way in.
  VectorXd theta_merr = theta_merr_in;
  if (int mode = merr_nig_mode())
    theta_merr.head(3) = nig_std::to_native(mode, theta_merr_in.head(3));
  if (!fix_flag[block_fix_theta_mu])
    theta_mu = theta_merr.segment(0, n_theta_mu);
  if (!fix_flag[block_fix_theta_sigma]) {
    int pos = 0;
    for (int i = 0; i < theta_sigma.size(); ++i) {
      if (!fix_theta_sigma_vec[i]) {
        theta_sigma(i) = theta_merr(n_theta_mu + pos);
        pos += 1;
      }
    }
  }
  if (!fix_flag[block_fix_theta_nu])
    theta_nu = theta_merr.segment(n_theta_mu + n_theta_sigma, n_theta_nu);

  // Cap nu at 1e4 (theta_nu is log nu).
  if (family != "normal" && theta_nu(0) > log(1e4))
    theta_nu(0) = log(1e4);

  // update rho, and Q_eps
  if (corr_measure && !fix_flag[block_fix_rho]) {
    rho(0) = th2rho(theta_merr(n_merr - 1));
    update_Q_eps(rho(0));
  }

  // update mu, sigma
  noise_mu = (B_mu * theta_mu);
  noise_sigma = (B_sigma * theta_sigma).array().exp();
  noise_nu = nu_lower_bound + (B_nu * theta_nu).array().exp();

  // show the Q construction
  // std::cout << "Q_eps == \n" << Q_eps << std::endl;
  // std::cout << "dQ_eps == \n" << dQ_eps << std::endl;
}

void BlockModel::sample_cond_noise_V(bool posterior) {
  if (family == "normal" || fix_flag[blcok_fix_V])
    return;
  // noise_V enters QQ and H through the measurement precision, so redrawing it
  // makes them stale. (Gaussian families return above and keep them.)
  invalidate_measurement();
  noise_prevV = noise_V;

  if (posterior) {
    VectorXd a_inc_vec = noise_mu.cwiseQuotient(noise_sigma).array().pow(2);
    VectorXd b_inc_vec = (get_residual() + noise_V.cwiseProduct(noise_mu))
                             .cwiseQuotient(noise_sigma)
                             .array()
                             .pow(2);
    // VectorXd b_inc_vec = (Y - A * getW() - X * beta -
    // noise_mu).cwiseQuotient(noise_sigma).array().pow(2);
    VectorXd a_vec_new = a_vec + a_inc_vec;
    VectorXd b_vec_new = b_vec + b_inc_vec;

    if (!corr_measure) {
      double dim = 1;
      VectorXd p_vec_new = p_vec - VectorXd::Constant(n_obs, 0.5 * dim);
      // noise_V = rGIG_cpp(p_vec_new, a_vec_new, b_vec_new, rng());
      NoiseUtil::sample_V(noise_V, family, p_vec_new, a_vec_new, b_vec_new,
                          rng);
    } else {
      // with pmu and psigma
      // pmat * res ~ N(-mu + mu V, Q^-1 = M^-1 diag(V) M^-T)
      // assert(noise_type == "nig");
      int dim = 2;
      VectorXd p_vec_new = p_vec - VectorXd::Constant(n_obs, 0.5 * dim);
      // loop over 1..n, sample V_i
      int i = 0;
      while (i < n_obs) {
        // std::cout << " i = " << i << std::endl;
        // std::cout << " has_cor[i] = " << has_correlation[i] << std::endl;
        if (has_correlation[i]) {
          // means obs_i and obs_i+1 share the same V
          noise_V[i] = noise_V[i + 1] =
              rGIG_cpp(p_vec[i], a_vec[i], b_vec[i], rng());
          i += 2;
        } else {
          noise_V[i] = noise_V[i + 1] =
              rGIG_cpp(p_vec[i], a_vec[i], b_vec[i], rng());
          i++;
        }
      }
    }
  } else {
    // noise_V = rGIG_cpp(p_vec, a_vec, b_vec, rng());
    NoiseUtil::sample_V(noise_V, family, p_vec, a_vec, b_vec, rng);
  }
}

// Draw W from its prior given V
void BlockModel::sampleW_V() {
  if (n_latent == 0)
    return;
  for (int li = 0; li < n_latent; ++li) {
    const VectorXd mean_li = latents[li]->getMean();   // mu (V - h)
    const VectorXd sv_li = latents[li]->getSV();       // sigma^2 V
    VectorXd z = NoiseUtil::rnorm_vec(mean_li.size(), 0, 1, rng());
    VectorXd rhs = mean_li + sv_li.cwiseSqrt().cwiseProduct(z);

    const SparseMatrix<double, 0, int> &Kl = latents[li]->getK();
    Eigen::SparseLU<SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int>> lu;
    lu.analyzePattern(Kl);
    lu.factorize(Kl);
    if (lu.info() != Eigen::Success)
      Rcpp::stop("sampleW_V(): factorization of K failed");

    latents[li]->setW(lu.solve(rhs));
    // The prior analogue of E(W | Y, V) is E(W | V) = K^-1 mu (V - h);
    // leaving cond_W at its stored posterior value would leak the fit.
    latents[li]->set_cond_W(lu.solve(mean_li));
  }
}

// posterior
Rcpp::List BlockModel::sampling(int n, int n_burnin, bool posterior,
                                const SparseMatrix<double> &A) {
  std::vector<VectorXd> AWs;     // blockA * blockW
  std::vector<VectorXd> Ws;      // return ZW (concatenated across latents)
  std::vector<VectorXd> cond_Ws; // return Z * E(W|Y,V)
  std::vector<VectorXd> Vs;      // blockV
  std::vector<VectorXd> mn_Vs;   // measurement nosie V
  rao_blackwell = true;

  // burn_in() runs posterior sweeps, so it is only meaningful when the
  // draws that follow are posterior draws.
  if (posterior && !all_gaussian)
    burn_in(n_burnin);

  for (int i = 0; i < n; i++) {
    if (posterior) {
      if (!all_gaussian)
        sample_cond_V();
      sampleW_VY();
      sample_cond_noise_V(true);
    } else {
      sample_uncond_V();
      sampleW_V();
      sample_cond_noise_V(false);
    }

    // Collect A Z W for output
    if (n_latent > 0) {
      VectorXd AZW = VectorXd::Zero(n_obs);
      VectorXd ZW = VectorXd::Zero(W_sizes);
      VectorXd ZcW = VectorXd::Zero(W_sizes);
      int woff = 0;
      for (int li = 0; li < n_latent; ++li) {
        const auto &Ai = latents[li]->getA();
        const auto &Zi = latents[li]->getZ();
        const auto &Wi = latents[li]->getW();
        const auto &cWi = latents[li]->get_cond_W();
        VectorXd ZiWi = Zi * Wi;
        VectorXd Zi_cWi = Zi * cWi;
        AZW.noalias() += Ai * ZiWi;
        ZW.segment(woff, ZiWi.size()) = ZiWi;
        ZcW.segment(woff, Zi_cWi.size()) = Zi_cWi;
        woff += ZiWi.size();
      }
      AWs.push_back(AZW);
      Ws.push_back(ZW);
      cond_Ws.push_back(ZcW);
    } else {
      AWs.push_back(VectorXd::Zero(n_obs));
      Ws.push_back(VectorXd::Zero(W_sizes));
      cond_Ws.push_back(VectorXd::Zero(W_sizes));
    }
    Vs.push_back(getV());
    // mn_Vs.push_back(var.getV());
  }

  return Rcpp::List::create(Rcpp::Named("AW") = AWs, Rcpp::Named("W") = Ws,
                            Rcpp::Named("V") = Vs,
                            Rcpp::Named("cond_W") = cond_Ws);
}

// fix parameter if converge
// void BlockModel::check_converge(vector<bool>& converge) {
//   int pos = 0;
//   for (std::vector<std::shared_ptr<Latent>>::iterator it = latents.begin();
//   it != latents.end(); it++) {
//     int theta_len = (*it)->get_n_params();
//     vector<bool> sub_converge (converge.begin(), converge.begin() +
//     theta_len);
//     (*it)->check_converge(sub_converge);
//     pos += theta_len;
//   }
// }

MatrixXd BlockModel::get_preconditioner() {
  if (!last_precond_valid) {
    throw std::runtime_error(
        "last_precond_valid is false, return identity matrix");
  }
  return last_precond;
}

// Main function for computing
void BlockModel::compute_grad_and_hessian(bool with_precond, double eps) {
  if (debug)
    ngme_io::out() << "Start compute_grad_and_hessian" << std::endl;
  ngme_timing::Scope _gt(ngme_timing::grad_total_us());
  auto t_total_start = std::chrono::steady_clock::now();
  long long t_sampleV_ms = 0, t_sampleW_ms = 0, t_rbtrace_ms = 0;
  long long t_build_s_ms = 0, t_set_s_ms = 0, t_dZ_ms = 0, t_grad_ms = 0;
  long long t_prec_latent_ms = 0, t_prec_ZGN_ms = 0, t_prec_merr_ms = 0;

  // Apply any parked operator-side probe budget before anything reads a probe
  // block, so no cache is left sized for the previous budget.
  if (pending_k_budget_ > 0) {
    for (auto &l : latents)
      l->set_n_trace_iter(pending_k_budget_);
    pending_k_budget_ = -1;
  }

  bool do_precond = with_precond;
  // Use the eps provided by the caller to keep caches consistent across layers

  VectorXd latent_grad = VectorXd::Zero(n_la_params);
  VectorXd noise_grad = VectorXd::Zero(n_params - n_la_params);
  MatrixXd precond_sum = MatrixXd::Zero(n_params, n_params);
  int precond_count = 0;

  // RB: one pass using conditional W; Gibbs: n_gibbs passes using sampled W
  // Only collapse to a single RB pass for all-Gaussian models.
  // For non-Gaussian we still need Gibbs over V even if RB is enabled.
  bool rb_all_gauss = (all_gaussian && rao_blackwell);
  int n_pass = rb_all_gauss ? 1 : n_gibbs;
  ngme_counters::add(ngme_counters::gibbs_passes, n_pass);
  // Use conditional mean of W for building observation score when RB is
  // requested
  bool use_condW = rao_blackwell;

  if (rao_blackwell)
    assemble_dK();

  // For covariance of gradient samples (only meaningful when n_pass>1)
  MatrixXd grad_samples(n_params, std::max(1, n_pass));

  for (int i = 0; i < n_pass; ++i) {
    // Gibbs sampler
    if (rb_all_gauss) {
      if (i == 0) {
        auto t_sw = std::chrono::steady_clock::now();
        sampleW_VY(); // compute QQ using current Z and V
        t_sampleW_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - t_sw)
                            .count();
        auto t_rb = std::chrono::steady_clock::now();
        compute_rb_trace();
        adapt_trace_probes(); // RB trace once (depends on QQ)
        t_rbtrace_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - t_rb)
                            .count();
      }
    } else {
      auto t_sv = std::chrono::steady_clock::now();
      { ngme_timing::Scope _pv(ngme_timing::grad_V_us());
      // Avoid duplicating QQ factorization here; sampleW_VY() updates QQ
      sample_cond_V();
      // Both V blocks are drawn before W.
      sample_cond_noise_V(); }
      t_sampleV_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t_sv)
                          .count();
      auto t_sw = std::chrono::steady_clock::now();
      sampleW_VY(false);
      t_sampleW_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t_sw)
                          .count();
      if (rao_blackwell) {
        auto t_rb = std::chrono::steady_clock::now();
        compute_rb_trace();
        adapt_trace_probes();
        t_rbtrace_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - t_rb)
                            .count();
      }
    }

    // Build observation score s = A^T D r for Z-chain
    auto t_bs = std::chrono::steady_clock::now();
    VectorXd s_full;
    { ngme_timing::Scope _ps(ngme_timing::grad_score_us());
    if (!corr_measure) {
      VectorXd residual = get_residual(use_condW);
      VectorXd inv_noise_SV =
          noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
      // Parenthesised for the same reason as M in sampleW_VY: left-to-right
      // grouping would build a scaled transpose of the sparse A (n_obs x W)
      // before touching the vector. Scaling the vector first leaves one
      // sparse-transpose-times-vector.
      s_full = A.transpose() * (inv_noise_SV.asDiagonal() * residual);
    } else {
      VectorXd residual = get_residual(use_condW);
      // Likewise: (A^T * Q_eps) is a sparse-sparse product, where
      // A^T * (Q_eps * residual) is a matrix-vector product twice.
      s_full = A.transpose() * (Q_eps * residual);
    }
    }
    t_build_s_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t_bs)
                        .count();

    // Aggregate gradients (latent + measurement Z-chain)
    auto t_g = std::chrono::steady_clock::now();
    ngme_timing::Scope _pg(ngme_timing::grad_assemble_us());
    VectorXd current_grad = VectorXd::Zero(n_params);
    int pos = 0;   // parameter offset
    int woff2 = 0; // W-slice offset for s_full
    for (int li = 0; li < n_latent; ++li) {
      auto &L = latents[li];
      L->compute_grad_and_hessian(rao_blackwell, do_precond);
      int theta_len = L->get_n_params();
      VectorXd gi = L->get_grad(); // excludes Z-chain measurement term; RB
                                   // compensated below
      // RB compensation at Block level
      if (rao_blackwell) {
        int n_k = L->get_n_theta_K();
        if (n_k > 0 && li < (int)rb_trace_K_latent.size() &&
            rb_trace_K_latent[li].size() == n_k) {
          gi.head(n_k) += rb_trace_K_latent[li];
        }
        int n_mu = L->get_n_theta_mu();
        int n_sig = L->get_n_theta_sigma();
        if (n_sig > 0 && li < (int)rb_trace_sigma_latent.size() &&
            rb_trace_sigma_latent[li].size() == n_sig) {
          gi.segment(n_k + n_mu, n_sig) += rb_trace_sigma_latent[li];
        }
      }
      // Add measurement Z-chain gradient: g_Z(j) = (dZ_j W_i)^T s_i
      int n_k = L->get_n_theta_K();
      if (n_k > 0) {
        VectorXd Wi_loc = use_condW ? L->get_cond_W() : L->getW();
        VectorXd s_i = s_full.segment(woff2, Wi_loc.size());
        for (int j = 0; j < n_k; ++j) {
          const auto &dZ_j = L->get_dZ(j);
          if (dZ_j.rows() == Wi_loc.size() && dZ_j.cols() == Wi_loc.size() &&
              dZ_j.nonZeros() > 0) {
            VectorXd gvec = dZ_j * Wi_loc;
            double gZ = gvec.dot(s_i);
            gi(j) += gZ;
          }
        }
      }
      // Standardised NIG coordinates: the optimiser works in
      // t = (log sigma_marg, zeta, log eta). Everything above is native, so the
      // chain rule is applied here, once, on the complete native gradient --
      // the RB and dZ terms added just above are native too and must be
      // included. grad_t = J^T grad_native, exactly (pure change of variables).
      {
        MatrixXd Jstd = L->get_nig_std_jacobian();
        if (Jstd.size() > 0) {
          int off = L->get_n_theta_K();
          gi.segment(off, 3) = Jstd.transpose() * gi.segment(off, 3).eval();
        }
      }
      current_grad.segment(pos, theta_len) = gi;
      latent_grad.segment(pos, theta_len) += gi;
      pos += theta_len;
      woff2 += L->get_W_size();
    }
    VectorXd noise_g = grad_theta_merr();
    // Standardised NIG coordinates for the measurement noise: chain rule on
    // the complete native gradient, grad_t = J^T grad_native.
    if (int mode = merr_nig_mode()) {
      Eigen::Vector3d native(theta_mu(0), theta_sigma(0), theta_nu(0));
      noise_g.head(3) =
          nig_std::jacobian(mode, native).transpose() * noise_g.head(3).eval();
    }
    current_grad.segment(n_la_params, n_merr) = noise_g;
    noise_grad.head(n_merr) += noise_g;
    if (!fix_flag[block_fix_beta]) {
      VectorXd beta_g = grad_beta();
      current_grad.segment(n_la_params + n_merr, n_feff) = beta_g;
      noise_grad.tail(n_feff) += beta_g;
    }
    _pg.stop();
    t_grad_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - t_g)
                     .count();

    // Running (Welford) variance of each gradient component. Combined with the
    // Hutchinson probe variance below this separates the two noise sources:
    //   Var(total) = Var(Gibbs) + Var(Hutchinson),  Var(Hutch) = probe_var / N.
    if (grad_diff_sq_.size() != n_params) {
      grad_prev_ = current_grad;
      grad_diff_sq_ = VectorXd::Zero(n_params);
      grad_run_n_ = 0;
    } else {
      // Exponentially weighted, never reset: a fresh window after every change
      // is short and noisy, which is what made the budget oscillate.
      const double a = 0.005; // ~200-iteration memory
      VectorXd d = current_grad - grad_prev_;
      if (grad_run_n_ == 0)
        grad_diff_sq_ = d.cwiseProduct(d);
      else
        grad_diff_sq_ = (1.0 - a) * grad_diff_sq_ + a * d.cwiseProduct(d);
      grad_prev_ = current_grad;
      grad_run_n_ += 1;
    }

    // Store gradient sample for covariance
    if (n_pass > 1)
      grad_samples.col(i) = current_grad;

    // Aggregate preconditioner blocks if requested
    if (do_precond) {
      // per-latent block preconditioners
      auto t_pl = std::chrono::steady_clock::now();
      ngme_timing::Scope _pl2(ngme_timing::grad_prec_lat_us());
      int pos2 = 0;
      // Measurement theta_sigma: the marginal Fisher information
      // (fisher_theta_sigma) or the complete-data Hessian, per
      // control_opt(precond_meas_sigma). "auto" uses Fisher for non-Gaussian
      // measurement noise only: from reasonable starting values the
      // complete-data block converges as fast for Gaussian noise at no extra
      // cost, while for non-Gaussian noise it stalls or diverges -- so
      // "complete" falls back to Fisher there. Fisher needs uncorrelated noise.
      // When it is used the complete-data cross terms are dropped: those with
      // beta and mu are exactly zero under the Gaussian Fisher information, and
      // keeping any beside a near-zero F_sigma could make the preconditioner
      // indefinite.
      const bool fisher_wanted = precond_meas_sigma_ == 1 || family != "normal";
      const bool fisher_sigma = fisher_wanted && !corr_measure &&
                                n_theta_sigma > 0 &&
                                !fix_flag[block_fix_theta_sigma];
      bool need_beta_cross = (n_feff > 0 && !fix_flag[block_fix_beta]);
      VectorXd inv_noise_SV_beta;
      if (need_beta_cross) {
        inv_noise_SV_beta =
            noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
      }

      for (int li = 0; li < n_latent; ++li) {
        int theta_len = latents[li]->get_n_params();
        MatrixXd Pi;
        try {
          Pi = latents[li]->preconditioner();
        } catch (...) {
          Pi = VectorXd::Constant(theta_len, latents[li]->get_V_size())
                   .asDiagonal();
        }
        precond_sum.block(pos2, pos2, theta_len, theta_len) += Pi;
        pos2 += theta_len;
      }
      _pl2.stop();
      t_prec_latent_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                              std::chrono::steady_clock::now() - t_pl)
                              .count();

      // Z-chain Hessian (measurement part):
      // H_{jk} = (Z_{jk} W)^T A^T D e  - (dZ_k W)^T A^T D A (dZ_j W)
      // We add the full (j,k) matrix per latent. If Z_{jk} is unavailable, we
      // fall back to the Gauss–Newton term (second term only).
      auto t_pz = std::chrono::steady_clock::now();
      ngme_timing::Scope _pz2(ngme_timing::grad_prec_ZGN_us());
      for (int li = 0; li < n_latent; ++li) {
        SparseMatrix<double> ADA_i;
        const auto &Ai = latents[li]->getA();
        if (!corr_measure) {
          VectorXd inv_noise_SV =
              noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
          ADA_i = Ai.transpose() * inv_noise_SV.asDiagonal() * Ai;
        } else {
          if (Q_eps.rows() == 0) {
            VectorXd inv_noise_SV =
                noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
            ADA_i = Ai.transpose() * inv_noise_SV.asDiagonal() * Ai;
          } else {
            ADA_i = Ai.transpose() * Q_eps * Ai;
          }
        }
        // Compute parameter and W offsets up to latent li
        int pbase = 0, woff_local = 0;
        for (int k = 0; k < li; ++k) {
          pbase += latents[k]->get_n_params();
          woff_local += latents[k]->get_W_size();
        }
        VectorXd Wi_loc =
            use_condW ? latents[li]->get_cond_W() : latents[li]->getW();
        VectorXd s_i = s_full.segment(woff_local, Wi_loc.size());
        int n_k = latents[li]->get_n_theta_K();
        if (n_k > 0) {
          MatrixXd HZ = MatrixXd::Zero(n_k, n_k);
          // Precompute dZ_j W and A dZ_j W vectors
          std::vector<VectorXd> dZW(n_k);
          std::vector<VectorXd> AZdZW(n_k);
          for (int j = 0; j < n_k; ++j) {
            const auto &dZ_j = latents[li]->get_dZ(j);
            if (dZ_j.rows() == Wi_loc.size() && dZ_j.cols() == Wi_loc.size() &&
                dZ_j.nonZeros() > 0) {
              dZW[j] = dZ_j * Wi_loc;
              AZdZW[j] = Ai * dZW[j];
            } else {
              dZW[j] = VectorXd::Zero(Wi_loc.size());
              AZdZW[j] = VectorXd::Zero(Ai.rows());
            }
          }
          // Build HZ
          for (int j = 0; j < n_k; ++j) {
            for (int k = 0; k < n_k; ++k) {
              // Second-derivative term: (Z_{jk} W)^T s_i ;
              double t1 = 0.0;
              const auto &d2Z_jk = latents[li]->get_d2Z(j, k);
              if (d2Z_jk.rows() == Wi_loc.size())
                t1 = (d2Z_jk * Wi_loc).dot(s_i);
              // GN term: − (dZ_k W)^T A^T D A (dZ_j W)
              double t2 = 0.0;
              if (dZW[k].size() > 0 && dZW[j].size() > 0)
                t2 = dZW[k].dot(ADA_i * dZW[j]);
              HZ(j, k) += (t1 - t2);
            }
          }
          // Accumulate Z Hessian into global preconditioner
          if (pbase + n_k <= n_params) {
            precond_sum.block(pbase, pbase, n_k, n_k) += HZ;
          }

          // (No complete-data cross term with measurement theta_sigma: see
          // fisher_sigma above. Its old formula also assumed diagonal Sigma.)

          // Cross-term with theta_mu (noise mean):
          // H_{theta,mu} = - J_theta^T diag(w_mu) B_mu,
          //   where w_mu = (V'-1)/(V' sigma^2) for uncorrelated
          if (n_theta_mu > 0 && !fix_flag[block_fix_theta_mu]) {
            int mu_col0 = n_la_params; // mu block starts the measurement part
            // Build w_mu per observation
            VectorXd w_mu =
                (noise_V.array() - 1.0).matrix().cwiseQuotient(noise_V);
            // include 1/sigma^2 factor (uncorrelated case uses diag(1/(sigma^2
            // V)))
            w_mu = w_mu.cwiseQuotient(noise_sigma.array().square().matrix());
            for (int j = 0; j < n_k; ++j) {
              const VectorXd &AZWj = AZdZW[j];
              if (AZWj.size() == n_obs) {
                VectorXd w = w_mu.cwiseProduct(AZWj);
                VectorXd row =
                    -(B_mu.transpose() * w); // 1 x n_theta_mu (as column)
                if (pbase + j < n_params && mu_col0 + n_theta_mu <= n_params) {
                  precond_sum.block(pbase + j, mu_col0, 1, n_theta_mu) +=
                      row.transpose();
                  precond_sum.block(mu_col0, pbase + j, n_theta_mu, 1) += row;
                }
              }
            }
          }

          // Cross-term with beta: H_{theta,beta} = - J_theta^T Sigma^{-1} X
          if (need_beta_cross) {
            int beta_col0 = n_la_params + n_merr;
            for (int j = 0; j < n_k; ++j) {
              const VectorXd &AZWj = AZdZW[j];
              if (AZWj.size() == n_obs) {
                VectorXd weighted = inv_noise_SV_beta.cwiseProduct(AZWj);
                VectorXd row = -(X.transpose() * weighted);
                if (pbase + j < n_params && beta_col0 + n_feff <= n_params) {
                  precond_sum.block(pbase + j, beta_col0, 1, n_feff) +=
                      row.transpose();
                  precond_sum.block(beta_col0, pbase + j, n_feff, 1) += row;
                }
              }
            }
          }
        }
      }
      _pz2.stop();
      t_prec_ZGN_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                           std::chrono::steady_clock::now() - t_pz)
                           .count();

      // measurement/fixed-effects block
      auto t_pm = std::chrono::steady_clock::now();
      ngme_timing::Scope _pm2(ngme_timing::grad_prec_merr_us());
      // (Measurement mu is a mean parameter: its block, and its cross term with
      // beta, are built jointly with beta in 1d below.)

      // 1b) Measurement sigma. Uncorrelated noise: the marginal Fisher
      // information (fisher_theta_sigma), with no mu-sigma cross term.
      if (fisher_sigma) {
        // It changes smoothly and each refresh costs a block of QQ solves, so it
        // is built at most once per iteration (in the first Gibbs pass) and
        // refreshed every fisher_refresh_every_ iterations, or sooner once
        // theta_sigma has moved by more than 0.05 -- it scales roughly like
        // sigma^4, so a stale value matters most while sigma is moving.
        bool refresh = fisher_sigma_cache_.rows() != n_theta_sigma ||
                       fisher_cache_theta_.size() != theta_sigma.size();
        if (!refresh && i == 0) {
          refresh = curr_iter - fisher_cache_iter_ >= fisher_refresh_every_ ||
                    (theta_sigma - fisher_cache_theta_).cwiseAbs().maxCoeff() >
                        0.05;
        }
        if (refresh) {
          fisher_sigma_cache_ = fisher_theta_sigma();
          fisher_cache_iter_ = curr_iter;
          fisher_cache_theta_ = theta_sigma;
        }
        precond_sum.block(n_la_params + n_theta_mu, n_la_params + n_theta_mu,
                          n_theta_sigma, n_theta_sigma) -= fisher_sigma_cache_;
      } else if (n_theta_sigma > 0 && !fix_flag[block_fix_theta_sigma]) {
        // Complete-data Hessian and mu-sigma cross: Gaussian measurement noise
        // under precond_meas_sigma "auto" or "complete", and correlated noise.
        // e = residual = Y - mu' (V-1) - A Z W - X beta
        VectorXd e = get_residual(use_condW);
        // H_sigma = -2 B_sigma^T diag(e^2 / (sigma^2 V)) B_sigma
        VectorXd wsig =
            2.0 *
            e.array().square().matrix().cwiseQuotient(
                noise_sigma.array().square().matrix().cwiseProduct(noise_V));
        MatrixXd Hsigma_full =
            -(B_sigma.transpose() * wsig.asDiagonal() * B_sigma);
        const std::vector<int> sig_free =
            free_sigma_cols(fix_theta_sigma_vec, (int)B_sigma.cols());
        // Place H_sigma just after the mu block within measurement corner
        for (int a = 0; a < (int)sig_free.size() && a < n_theta_sigma; ++a)
          for (int b = 0; b < (int)sig_free.size() && b < n_theta_sigma; ++b)
            precond_sum(n_la_params + n_theta_mu + a,
                        n_la_params + n_theta_mu + b) +=
                Hsigma_full(sig_free[a], sig_free[b]);
        // Cross H_{mu,sigma} = -2 B_mu^T diag(((V-1) ⊙ e) / (sigma^2 V))
        // B_sigma
        if (n_theta_mu > 0 && !fix_flag[block_fix_theta_mu]) {
          VectorXd wms =
              2.0 * (noise_V - VectorXd::Ones(n_obs))
                        .cwiseProduct(e)
                        .cwiseQuotient(
                            noise_sigma.array().square().matrix().cwiseProduct(
                                noise_V));
          MatrixXd Hmu_sigma_full =
              -(B_mu.transpose() * wms.asDiagonal() * B_sigma);
          for (int a = 0; a < (int)sig_free.size() && a < n_theta_sigma; ++a)
            for (int m = 0; m < n_theta_mu; ++m) {
              const double v = Hmu_sigma_full(m, sig_free[a]);
              precond_sum(n_la_params + m, n_la_params + n_theta_mu + a) += v;
              precond_sum(n_la_params + n_theta_mu + a, n_la_params + m) += v;
            }
        }
      }

      // 1c) Analytic Hessian for measurement nu (no cross-terms)
      if (n_theta_nu > 0 && !fix_flag[block_fix_theta_nu]) {
        if (family != "normal") {
          // Use the same analytic form as latent: H_nu = - B_nu^T diag(nu ⊙ c)
          // B_nu for NIG, and the appropriate GAL/t variants handled inside
          // NoiseUtil.
          MatrixXd Hnu = -NoiseUtil::hess_theta_nu(
              family, B_nu, noise_nu, noise_V, VectorXd::Ones(noise_V.size()),
              nu_lower_bound);
          precond_sum.block(n_la_params + n_theta_mu + n_theta_sigma,
                            n_la_params + n_theta_mu + n_theta_sigma,
                            n_theta_nu, n_theta_nu) += Hnu;
        }
      }

      // H_{sigma,beta} = -2 B_sigma^T diag(e/(sigma^2 ∘ V)) X, complete-data;
      // zero under the Fisher information, so only with the complete-data block.
      if (!fisher_sigma && n_theta_sigma > 0 && n_feff > 0 &&
          !fix_flag[block_fix_theta_sigma] && !fix_flag[block_fix_beta]) {
        VectorXd e_sb = get_residual(use_condW);
        VectorXd wsb =
            2.0 *
            e_sb.cwiseQuotient(
                noise_sigma.array().square().matrix().cwiseProduct(noise_V));
        MatrixXd Hsigma_beta_full =
            -(B_sigma.transpose() * wsb.asDiagonal() * X);
        const std::vector<int> sb_free =
            free_sigma_cols(fix_theta_sigma_vec, (int)B_sigma.cols());
        for (int a = 0; a < (int)sb_free.size() && a < n_theta_sigma; ++a)
          for (int b = 0; b < n_feff; ++b) {
            const double v = Hsigma_beta_full(sb_free[a], b);
            precond_sum(n_la_params + n_theta_mu + a, n_la_params + n_merr + b)
                += v;
            precond_sum(n_la_params + n_merr + b, n_la_params + n_theta_mu + a)
                += v;
          }
      }

      // 1d) Mean parameters -- fixed effects beta and measurement mu -- jointly,
      // Rao-Blackwellised over W. Given V both enter only the mean of Y,
      //   E[Y | W, V] = diag(V - 1) B_mu theta_mu + X beta + A Z W,
      // and their scores are linear in W, so with M = [diag(V - 1) B_mu, X]
      // Louis' identity gives, exactly given V,
      //   H = -(M^T D M - B^T QQ^{-1} B),   B = (AZ)^T D M,
      // D the measurement precision (Q_eps when correlated): in the Gaussian case
      // the marginal -M^T Sigma_Y^{-1} M. The complete-data -M^T D M overstates
      // it by the information W carries -- badly for an intercept or a smooth
      // covariate under a correlated field, and without bound for mu, whose
      // sum (V-1)^2 / (sigma^2 V) diverges as sigma -> 0 while the marginal
      // stays finite. Mean parameters are Fisher-orthogonal to sigma, so there
      // are no mu-sigma or beta-sigma cross terms. One QQ solve per parameter.
      const bool mean_mu = n_theta_mu > 0 && !fix_flag[block_fix_theta_mu];
      const bool mean_beta = n_feff > 0 && !fix_flag[block_fix_beta];
      if (mean_mu || mean_beta) {
        const int p_mu = mean_mu ? n_theta_mu : 0;
        const int p_b = mean_beta ? n_feff : 0;
        MatrixXd M(n_obs, p_mu + p_b);
        if (mean_mu)
          M.leftCols(p_mu) =
              (noise_V.array() - 1.0).matrix().asDiagonal() * B_mu;
        if (mean_beta)
          M.rightCols(p_b) = X;
        MatrixXd DM;
        if (!corr_measure)
          DM = meas_prec().asDiagonal() * M;
        else
          DM = Q_eps * M;
        MatrixXd Hm = -(M.transpose() * DM);
        if (n_latent > 0 && QQ_valid) {
          MatrixXd B = get_AZ().transpose() * DM;
          MatrixXd QQinvB(B.rows(), B.cols());
          for (int l = 0; l < B.cols(); ++l) {
            VectorXd b = B.col(l);
            QQinvB.col(l) = chol_QQ.solve(b);
          }
          Hm += B.transpose() * QQinvB;
          Hm = (0.5 * (Hm + Hm.transpose())).eval();
        }
        // Prior curvature, by central difference of the prior score.
        auto prior_curv = [](const string &type, const VectorXd &param,
                             double v) {
          const double h = 1e-5 * std::max(1.0, std::abs(v));
          return (PriorUtil::d_log_dens(type, param, v + h) -
                  PriorUtil::d_log_dens(type, param, v - h)) /
                 (2.0 * h);
        };
        if (mean_mu && prior_mu_target == "coef")
          for (int l = 0; l < p_mu; ++l)
            Hm(l, l) += prior_curv(prior_mu_type, prior_mu_param, theta_mu(l));
        for (int l = 0; l < p_b; ++l)
          Hm(p_mu + l, p_mu + l) +=
              prior_curv(prior_beta_type[l], prior_beta_param[l], beta(l));

        const int mu0 = n_la_params, b0 = n_la_params + n_merr;
        if (mean_mu)
          precond_sum.block(mu0, mu0, p_mu, p_mu) += Hm.topLeftCorner(p_mu, p_mu);
        if (mean_beta)
          precond_sum.block(b0, b0, p_b, p_b) += Hm.bottomRightCorner(p_b, p_b);
        if (mean_mu && mean_beta) {
          precond_sum.block(mu0, b0, p_mu, p_b) += Hm.topRightCorner(p_mu, p_b);
          precond_sum.block(b0, mu0, p_b, p_mu) +=
              Hm.bottomLeftCorner(p_b, p_mu);
        }
      }

      _pm2.stop();
      t_prec_merr_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - t_pm)
                            .count();

      precond_count += 1;
    }
    // std::cout << "precond_sum = " << precond_sum << std::endl;
  }

  // Gradient covariance across samples (only if >1)
  if (n_pass <= 1) {
    grad_covariance = MatrixXd::Zero(n_params, n_params);
  } else {
    VectorXd grad_mean = grad_samples.rowwise().mean();
    MatrixXd centered = grad_samples.colwise() - grad_mean;
    grad_covariance = (1.0 / (n_pass - 1)) * centered * centered.transpose();
  }

  // Average accumulated gradients
  latent_grad /= std::max(1, n_pass);
  noise_grad /= std::max(1, n_pass);

  VectorXd avg_gradient = VectorXd::Zero(n_params);
  avg_gradient.head(n_la_params) = latent_grad;
  avg_gradient.tail(n_params - n_la_params) = noise_grad;

  // Finalize cached preconditioner
  if (do_precond) {
    if (precond_count > 0) {
      last_precond = (1.0 / precond_count) * precond_sum;
      // Standardised measurement NIG coordinates: the blocks above are native,
      // so H_t = T^T H T with T the identity except J = d(native)/dt on the
      // measurement (mu, sigma, nu) block. Cross blocks pick up J on that side.
      if (int mode = merr_nig_mode()) {
        Eigen::Vector3d native(theta_mu(0), theta_sigma(0), theta_nu(0));
        MatrixXd T = MatrixXd::Identity(n_params, n_params);
        T.block(n_la_params, n_la_params, 3, 3) = nig_std::jacobian(mode, native);
        last_precond = (T.transpose() * last_precond * T).eval();
      }
      // last_precond is the Hessian; Ngme::precond() negates it into the
      // information. Subtract so the ridge actually increases the
      // information diagonal that llt() factorises.
      last_precond -= VectorXd::Constant(n_params, 1e-5).asDiagonal();
      last_precond_valid = true;
    } else {
      last_precond_valid = false;
    }
  }

  last_gradient = avg_gradient;
  last_grad_valid = true;

  // Publish the sub-phases so an iteration can be accounted whether or not
  // debug printing is on. Millisecond resolution is what the existing counters
  // carry; over a fit that is far finer than the numbers being compared.
  ngme_timing::add(ngme_timing::grad_sampleV_us(), t_sampleV_ms * 1000);
  ngme_timing::add(ngme_timing::grad_sampleW_us(), t_sampleW_ms * 1000);
  ngme_timing::add(ngme_timing::grad_rbtrace_us(), t_rbtrace_ms * 1000);

  if (debug) {
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t_total_start)
                        .count();
    ngme_io::out() << "[block] compute_grad_and_hessian timing (ms): total="
                << total_ms << ", sampleV=" << t_sampleV_ms
                << ", sampleW=" << t_sampleW_ms << ", rb_trace=" << t_rbtrace_ms
                << ", build_s=" << t_build_s_ms << ", set_s=" << t_set_s_ms
                << ", dZ=" << t_dZ_ms << ", grad=" << t_grad_ms
                << ", prec_latent=" << t_prec_latent_ms
                << ", prec_ZGN=" << t_prec_ZGN_ms
                << ", prec_merr=" << t_prec_merr_ms << std::endl;
  }
}

// Size the Hutchinson probe budget so the trace estimator's contribution to the
// gradient variance is a small fraction of the Gibbs sampling noise it sits on
// top of. Var(Hutch) = probe_var / N falls as 1/N, so the N that reaches
// Var(Hutch) = frac * Var(Gibbs) is  N = probe_var / (frac * Var(Gibbs)).
// Taking the max over the affected parameters keeps the worst one in budget.
// Probe variance is measured within an iteration and the total across
// iterations, so the two sources are separable without extra work.
void BlockModel::adapt_trace_probes() {
  if (!rao_blackwell || !trace_adapt || grad_run_n_ < 30)
    return;
  // Under the cost rule this function only keeps the statistics current. The
  // search runs at a fixed budget -- its job is to reach stationarity, and
  // probe noise was measured not to change how long that takes -- and the
  // budget is chosen by suggest_trace_N_cost() at the polish checkpoints. The
  // share arithmetic below would otherwise run to completion once the polish
  // set cost_ratio_, and write a suggestion nothing ever reads.
  if (trace_adapt_cost_rule) {
    update_probe_var_stats();
    return;
  }
  // The denominator scales the probe budget so that the
  // Hutchinson variance is trace_adapt_frac of the total per-iteration gradient
  // variance, but for an all-Gaussian model with Rao-Blackwell there is a
  // single pass using the conditional mean of W, so the probes are the ONLY
  // randomness and the share being targeted is identically one. What keeps the
  // budget finite is that var_tot is E[(g_t-g_{t-1})^2]/2, which is a variance
  // only while the gradient's mean is constant and otherwise also contains the
  // systematic drift of the gradient as the optimizer moves. So the rule really
  // drives probe noise below a fraction of the squared per-iteration drift, and
  // as a fit converges that drift goes to zero and the budget climbs to
  // trace_adapt_max. Whether more probes are worth their cost depends on the
  // operator, and for the tested cases, it seems like the probes are so cheap
  // that it is worth it as it reduces the total number of iterations needed.
  // No countdown here. This runs every pass so the smoothed statistics stay
  // current and the suggestion tracks them; WHEN a suggestion is acted on is
  // the driver's business, and it acts at the convergence checkpoints where the
  // chains can be given one budget between them. Keeping a countdown here as
  // well would put a second, invisible cadence behind that one -- and this one
  // counted CALLS, of which there is one per Gibbs pass, so it ran faster for a
  // non-Gaussian noise than for a Gaussian one without saying so.
  if (chol_QQ.is_exact_trace())
    return; // already exact, nothing to gain

  // Per-iteration gradient variance, robust to the mean drifting under the
  // optimizer: E[(g_t - g_{t-1})^2] / 2.
  VectorXd var_tot = grad_diff_sq_ / 2.0;

  update_probe_var_stats();
  const VectorXd &pv = probe_var_ewma_;

  const int N_cur = std::max(1, chol_QQ.get_N_iter());

  // Work with the BOUNDED share r = Var(Hutchinson) / Var(total) in [0,1)
  // rather than solving for N directly.
  double r_max = 0.0;
  for (int j = 0; j < n_params; ++j) {
    if (!(pv(j) > 0.0) || !(var_tot(j) > 0.0))
      continue;
    double r = (pv(j) / (double)N_cur) / var_tot(j);
    if (r > 0.95)
      r = 0.95; // beyond this the split is not measurable; do not extrapolate
    r_max = std::max(r_max, r);
  }
  if (!(r_max > 0.0))
    return;

  const double r_t = std::min(std::max(trace_adapt_frac, 1e-3), 0.9);
  // Deadband: leave the budget alone while the share is in the right region.
  // Without it the estimate's own noise is enough to keep moving N every time.
  if (r_max > 0.5 * r_t && r_max < 2.0 * r_t)
    return;

  // Var(Hutch) falls as 1/N, so the exact multiplier for r_max -> r_t is
  //   N_new/N_cur = [r_max/(1-r_max)] * [(1-r_t)/r_t].
  double factor = (r_max / (1.0 - r_max)) * ((1.0 - r_t) / r_t);
  // Take only a fraction of that step, in log space, and cap it. The budget
  // then approaches its target geometrically over several updates instead of
  // jumping to the solved value on one reading and bouncing off the opposite
  // bound next time.
  factor = std::pow(factor, 0.3);
  factor = std::min(std::max(factor, 0.8), 1.25);

  int N_new = (int)std::lround(N_cur * factor);
  N_new = std::min(std::max(N_new, trace_adapt_min), trace_adapt_max);
  // No floor is applied here. The solver clamps its own budget to the probing
  // floor (sparse_llt_solver::refresh_budget_), so a suggestion below it simply
  // has no effect and the boundary cannot be cycled across.
  // Suggest, do not apply: parallel chains have to end up on the SAME budget,
  // and the only place they are all stopped together is the convergence
  // checkpoint, so the driver decides there. See suggested_trace_N_.
  if (N_new != N_cur)
    suggested_trace_N_ = N_new;
}

// Gather the raw probe variances into the parameter layout and smooth them.
// Each reading is a sample variance over as few as five probes, far too noisy
// to steer on unsmoothed -- that was the main cause of budget oscillation. Both
// rules read the smoothed vector, and the cost rule needs it kept current
// through the search while not acting on it, hence the split.
void BlockModel::update_probe_var_stats() {
  VectorXd pv = VectorXd::Zero(n_params);
  int pos = 0;
  for (int li = 0; li < n_latent; ++li) {
    int n_k = latents[li]->get_n_theta_K();
    int n_mu = latents[li]->get_n_theta_mu();
    int n_sig = latents[li]->get_n_theta_sigma();
    if (li < (int)rb_probe_var_K_latent.size() &&
        rb_probe_var_K_latent[li].size() == n_k)
      pv.segment(pos, n_k) = rb_probe_var_K_latent[li];
    if (li < (int)rb_probe_var_sigma_latent.size() &&
        rb_probe_var_sigma_latent[li].size() == n_sig)
      pv.segment(pos + n_k + n_mu, n_sig) = rb_probe_var_sigma_latent[li];
    pos += latents[li]->get_n_params();
  }
  if (rb_probe_var_noise_sigma.size() == n_theta_sigma)
    pv.segment(n_la_params + n_theta_mu, n_theta_sigma) =
        rb_probe_var_noise_sigma;
  if (probe_var_ewma_.size() != n_params)
    probe_var_ewma_ = pv;
  else
    probe_var_ewma_ = 0.9 * probe_var_ewma_ + 0.1 * pv;
}

// Size the probe budget by what it costs, not by the share of variance it
// carries. Used from the polish on, where the reported estimate is the
// Polyak-Ruppert average of the iterates.
//
// For that average Var(theta_bar_T) goes as V(N)/T, and a pass costs
// c(N) = a + bN, so a fixed amount of work buys T = C/c(N) iterations and the
// precision goes as V(N)c(N)/C. Minimising
//     V(N) c(N) = (V_gibbs + P/N)(a + bN)
// gives N* = sqrt((P / V_gibbs) * (a / b)), with P the raw probe variance
// (Var(Hutchinson) = P/N) and V_gibbs the noise that is not from the probes.
// Both come from quantities the share rule already computes; a/b is the cost
// ratio taken once when the polish begins.
//
// The share rule instead targets Var(Hutch) = frac * var_tot, where var_tot
// itself contains Var(Hutch) -- so it reduces to holding probe noise to a
// fraction of the drift, which goes to zero on a converging fit, sending the
// budget to its cap whatever it measures. This rule has no such hole, and it
// responds to a better estimator by asking for fewer probes: structured probing
// cuts P several-fold, so N* falls by the square root of that. Where V_gibbs is
// negligible N* is unbounded in principle, but the objective saturates once
// bN >> a -- past there a probe and an iteration buy the same thing, which is
// what the cap below is.
int BlockModel::suggest_trace_N_cost() {
  if (!(cost_ratio_ > 0.0) || chol_QQ.is_exact_trace())
    return -1;
  if (probe_var_ewma_.size() != n_params || grad_diff_sq_.size() != n_params)
    return -1;
  const int N_cur = std::max(1, chol_QQ.get_N_iter());
  const VectorXd var_tot = grad_diff_sq_ / 2.0;

  // rho = P / V_gibbs, per parameter. A ratio of two variances of the same
  // parameter, so it is free of that parameter's scale and the values can be
  // compared across them without any standardisation.
  double rho_max = 0.0;
  for (int j = 0; j < n_params; ++j) {
    const double P = probe_var_ewma_(j);
    if (!(P > 0.0) || !R_finite(var_tot(j)) || !(var_tot(j) > 0.0))
      continue;
    const double v_probe = P / (double)N_cur;
    // What is left once the probes are accounted for. Floored well away from
    // zero: var_tot is itself an estimate, so the difference can come out
    // negative or minutely positive on a fit where the probes carry everything,
    // and rho would then be enormous for no measurable reason. The floor caps
    // rho at 1/eps and the saturation cap below takes over from there.
    const double v_gibbs = std::max(var_tot(j) - v_probe, 1e-3 * v_probe);
    rho_max = std::max(rho_max, v_probe * (double)N_cur / v_gibbs);
  }
  if (!(rho_max > 0.0))
    return -1;

  const double n_star = std::sqrt(rho_max * cost_ratio_);
  // Saturation: past the point where the probes cost as much as everything
  // else, another probe and another iteration buy the same thing.
  int N_new = (int)std::lround(std::min(n_star, cost_ratio_));
  N_new = std::min(std::max(N_new, trace_adapt_min), trace_adapt_max);
  // The polish only ever spends less: its job is to average, and another
  // iteration of averaging beats another probe.
  N_new = std::min(N_new, N_cur);
  // But never out of probing: rho was measured under the scheme running now, so
  // falling back to dense would give back the variance the colouring bought and
  // the chosen budget would not deliver what it was chosen for.
  if (chol_QQ.probing_active()) {
    const int floor_budget = chol_QQ.probe_min_budget();
    if (floor_budget > 0)
      N_new = std::max(N_new, floor_budget);
  }
  return N_new;
}

bool BlockModel::begin_polish_trace_rule() {
  if (!rao_blackwell || !trace_adapt || !trace_adapt_cost_rule)
    return false;

  // a/b -- probe columns per pass -- from the fill of QQ's factor, not a clock.
  // A probe column is two triangular solves, O(nnz(L)); the pass is dominated
  // by the factorization, O(sum_j |L(:,j)|^2); their ratio is of order nnz(L)/n.
  // That is a property of the matrix and its ordering, so the budget it implies
  // is the same in every run -- an elapsed-time ratio is not: the same search
  // twice gave 63 and 20, moving the polish between 60 probes and 10. The gate
  // for the selected inverse already computes the fill, so this is free.
  //
  // It underestimates a/b -- probe solves run nearer peak than the
  // factorization, and a pass holds sampling and assembly it does not cover --
  // so the rule asks for fewer probes than the optimum. The safe direction.
  if (!(qq_fill_ > 0.0))
    return false;
  cost_ratio_ = std::max(1.0, qq_fill_);

  // Steering needs the probe variance, which only whole replicates can measure
  // -- two draws over the colouring. Ask for that only if the budget affords
  // it: otherwise the search already sits at what a pass affords, and dropping
  // to dense to buy a number used only to hold the budget would give back the
  // variance the colouring was buying.
  if (chol_QQ.probing_active()) {
    const int p = chol_QQ.probe_colours();
    if (p < 1 || 2 * p > chol_QQ.get_N_iter()) {
      if (debug)
        ngme_io::out() << "[trace_cost] holding the search budget: two draws "
                          "over the colouring (" << 2 * p
                       << ") exceed what a pass affords ("
                       << (int)std::lround(cost_ratio_) << ")\n";
      return false;
    }
  }

  in_polish_ = true;
  setup_qq_probing();
  if (debug)
    ngme_io::out() << "[trace_cost] a pass costs about " << cost_ratio_
                   << " probe columns (fill); saturation budget "
                   << (int)std::lround(cost_ratio_) << "\n";
  return true;
}

void BlockModel::apply_trace_N(int N) {
  if (N < 1 || chol_QQ.is_exact_trace())
    return;
  suggested_trace_N_ = -1; // consumed
  if (N == chol_QQ.get_N_iter())
    return;
  if (debug)
    ngme_io::out() << "[trace_adapt] probes " << chol_QQ.get_N_iter() << " -> "
                   << N << "\n";
  chol_QQ.set_N_iter(N);
  // The operator-side budget has to move with it. cholK_solver is initialised
  // once, so a budget left behind here stays at its construction value for the
  // whole run however far the QQ budget travels -- and those probes feed the
  // gradient of theta_K directly. Parked rather than applied: the operator
  // caches a factored probe block across Gibbs passes, so it must not be
  // resized from outside that computation; update_all() picks the value up.
  if (trace_adapt_k)
    pending_k_budget_ = N;
  // The trace statistics are exponentially weighted and stay valid across the
  // change.
}

namespace {
// The solver reports a negative spread when it could not measure one -- probing
// with a single sign draw over the colouring. The budget controller treats a
// non-positive reading as "nothing to steer on" and leaves the budget alone,
// which is the right response, so map it onto zero here rather than letting a
// negative number into the smoothed statistics.
inline double probe_var_of(const sparse_llt_solver &s) {
  const double v = s.last_probe_var();
  return v > 0.0 ? v : 0.0;
}

// How many times the probe colouring may be re-sourced before probing is given
// up on. See setup_qq_probing().
constexpr int kMaxProbingSetups = 3;
} // namespace

void BlockModel::setup_qq_probing() {
  if (!trace_probing) {
    chol_QQ.disable_probing();
    return;
  }
  // The colouring is a per-pattern cost, amortized over every iteration that
  // shares the pattern -- which for almost every model is the whole fit, since
  // the symbolic phase runs once. A pattern that keeps moving never amortizes
  // it: a rational approximation rebuilds QQ's structure as its parameters
  // move, and re-colouring at each symbolic phase would cost about what the
  // probes it is trying to improve cost. Probing is dropped for the rest of
  // the fit rather than paid for on every iteration.
  if (++qq_probing_setups_ > kMaxProbingSetups) {
    chol_QQ.disable_probing();
    return;
  }
  // A spread can only be measured across whole replicates, so two sign draws
  // over the colouring are what it takes to report a probe variance -- and the
  // budget controller runs on that number. But only where a budget is actually
  // being steered: under the cost rule the search holds its budget, so it needs
  // no variance and one draw over one colouring will do. That is the cheapest
  // form of probing available, and it is what the search runs on.
  const bool need_probe_var = trace_adapt && (!trace_adapt_cost_rule || in_polish_);
  const int min_reps = need_probe_var ? 2 : 1;
  // No colouring larger than the largest budget this fit can reach is ever
  // usable, so the search is capped there rather than run to completion on a
  // graph too dense for this to pay.
  const int max_colours =
      std::max(chol_QQ.get_N_iter(), trace_adapt ? trace_adapt_max : 1);
  // Whether to pay the colouring floor. It is a hard floor: a colouring must be
  // used whole, and any partition coarser than distance-1 measures the same as
  // dense, so there is no cheaper form of probing to fall back on.
  //
  // For an all-Gaussian Rao-Blackwell fit the answer needs no measurement --
  // the gradient uses the conditional mean of W and nothing else is sampled, so
  // the probes are its only randomness and the floor is worth paying at once.
  // Elsewhere the Gibbs sampling carries noise the probes cannot remove, so the
  // caller's budget stands and the cost rule decides in the polish.
  const bool probes_are_the_only_noise = all_gaussian && rao_blackwell;
  const int requested = chol_QQ.get_requested_N_iter();
  const int raise_cap =
      trace_adapt_cost_rule
          ? ((!in_polish_ && probes_are_the_only_noise)
                 ? std::max(trace_adapt_max, requested)
                 : 0)
          : (trace_probing_raise_budget > 1.0
                 ? std::min((int)std::floor(trace_probing_raise_budget * requested),
                            trace_adapt_max)
                 : 0);
  chol_QQ.set_probing_source(QQ, trace_probing_max_dist, max_colours, min_reps,
                             raise_cap);
}

// Selected inversion is exact and, for a low-fill factor, cheaper than even a
// handful of probes. The fill ratio is the operational test and is decided once per fit.
double BlockModel::qq_trace(const SparseMatrix<double> &T, double &probe_var) {
  if (selinv_state_ < 0) {
    // Decided once per fit. Take the free lower bound on the fill first
    double fr = chol_QQ.fill_lower_bound();
    bool ruled_out = !chol_QQ.selinv_supported() || fr > selinv_max_fill;
    if (!ruled_out) {
      fr = chol_QQ.fill_ratio();
      ruled_out = fr > selinv_max_fill;
    }
    selinv_state_ = ruled_out ? 0 : 1;
    qq_fill_ = fr;
    if (debug)
      ngme_io::out() << "[selinv] fill_ratio=" << fr
                  << " -> " << (selinv_state_ ? "exact selected inverse"
                                            : "Hutchinson probes")
                << "\n";
    // Nothing on the selected-inverse path will be touched again, so hand back
    // everything it owns: its private factor, the symmetric copy of QQ it
    // factorizes, and the selected inverse itself.
    if (selinv_state_ == 0)
      chol_QQ.disable_selinv();
  }
  if (selinv_state_ == 1) {
    double v = 0.0;
    if (chol_QQ.selinv_trace(T, v)) {
      probe_var = 0.0; // exact: contributes no estimation variance
      return v;
    }
    selinv_state_ = 0; // pattern not covered; fall back for the rest of the fit
    chol_QQ.disable_selinv();
  }
  double v = chol_QQ.trace(T, rng());
  probe_var = probe_var_of(chol_QQ);
  report_probing();
  return v;
}

void BlockModel::report_probing() {
  if (!debug)
    return;
  const int d = chol_QQ.probing_active() ? chol_QQ.probe_dist() : 0;
  const int r = chol_QQ.probing_active() ? chol_QQ.probe_reps() : 0;
  if (d == probing_reported_dist_ && r == probing_reported_reps_)
    return;
  probing_reported_dist_ = d;
  probing_reported_reps_ = r;
  if (d > 0)
    ngme_io::out() << "[probing] distance=" << d
                   << " colours=" << chol_QQ.probe_colours() << " reps=" << r
                   << " probes=" << chol_QQ.get_n_probes() << " of budget "
                   << chol_QQ.get_N_iter() << "\n";
  else
    ngme_io::out() << "[probing] off: " << chol_QQ.get_n_probes()
                   << " dense Rademacher probes\n";
}

double BlockModel::qq_trace_factored_shared(const SparseMatrix<double> &A,
                                            const VectorXd &d,
                                            const SparseMatrix<double> &B,
                                            const Eigen::MatrixXd &BQU,
                                            bool have_BQU, double &probe_var) {
  if (selinv_state_ == 0 && have_BQU) {
    double v = chol_QQ.trace_factored_with(A, d, BQU);
    probe_var = probe_var_of(chol_QQ);
    report_probing();
    return v;
  }
  return qq_trace_factored(A, d, B, probe_var);
}

double BlockModel::qq_trace_factored(const SparseMatrix<double> &A,
                                     const VectorXd &d,
                                     const SparseMatrix<double> &B,
                                     double &probe_var) {
  // selinv_state_ is -1 until the fill test has run and 1 when the selected
  // inverse won it; both need the assembled matrix, so build it and take the
  // ordinary path. Once probing is settled on the product is never formed.
  static const bool force_assemble = [] {
    const char *e = std::getenv("NGME2_NO_FACTORED_TRACE");
    return e && *e && std::string(e) != "0" && std::string(e) != "false";
  }();
  if (selinv_state_ != 0 || force_assemble) {
    SparseMatrix<double> T = A.transpose() * d.asDiagonal() * B;
    return qq_trace(T, probe_var);
  }
  double v = chol_QQ.trace_factored(A, d, B, rng());
  probe_var = probe_var_of(chol_QQ);
  report_probing();
  return v;
}

void BlockModel::compute_rb_trace() {
  if (debug)
    ngme_io::out() << "start compute trace" << std::endl;
  auto t_start = std::chrono::steady_clock::now();
  int n = 0;
  int woff = 0; // offset in W-space for embedding latent-local pieces
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
  // Ensure storage is sized for all latents
  rb_trace_K_latent.resize(n_latent);
  rb_trace_sigma_latent.resize(n_latent);
  rb_probe_var_K_latent.resize(n_latent);
  rb_probe_var_sigma_latent.resize(n_latent);

  for (int i = 0; i < n_latent; i++) {
    VectorXd rb_trace_K(latents[i]->get_n_theta_K());
    VectorXd rb_trace_sigma(latents[i]->get_n_theta_sigma());
    VectorXd pv_K = VectorXd::Zero(latents[i]->get_n_theta_K());
    VectorXd pv_sigma = VectorXd::Zero(latents[i]->get_n_theta_sigma());


    // Both trace loops below use B = K, so K*QU is shared between them and
    // across every parameter instead of being rebuilt per call. Only on the
    // probe path: the selected inverse wants the assembled product.
    Eigen::MatrixXd KQU;
    bool have_KQU = false;
    if (selinv_state_ == 0) {
      KQU = chol_QQ.trace_factored_rhs(K, rng());
      have_KQU = true;
    }

    // compute for K: tr(QQ^-1 dK^T diag(1/SV) K)
    { ngme_timing::Scope _s(ngme_timing::rb_sec_K_us());
    for (int j = 0; j < latents[i]->get_n_theta_K(); j++) {
      { double pvj = 0.0;
        rb_trace_K[j] = -qq_trace_factored_shared(block_dK[i][j], inv_SV, K,
                                                  KQU, have_KQU, pvj);
        pv_K[j] += pvj; }
    }
    }

    // compute for sigma: tr(Q^-1 K B_sigma.col(j)/SV K^T) for non-fixed
    // theta_sigma
    ngme_timing::Scope _ss(ngme_timing::rb_sec_sigma_us());
    vector<bool> fix_theta_sigma_vec = latents[i]->get_theta_unfixed_sigma();
    // One free component per j, so consecutive fixed components all have to be
    // stepped over; skipping a single one lands on a fixed column whenever two
    // fixed components sit together.
    const std::vector<int> lat_sig_free =
        free_sigma_cols(fix_theta_sigma_vec, (int)fix_theta_sigma_vec.size());
    for (int j = 0; j < latents[i]->get_n_theta_sigma(); j++) {
      const int pos = (j < (int)lat_sig_free.size()) ? lat_sig_free[j] : j;

      // build B_sigma_col_j (consider all latents)
      VectorXd BSigma_col_over_SV = VectorXd::Zero(V_sizes);
      BSigma_col_over_SV.segment(n, latents[i]->get_V_size()) =
          latents[i]->get_BSigma_col(pos);
      BSigma_col_over_SV = BSigma_col_over_SV.cwiseProduct(inv_SV);

      { double pvj = 0.0;
        rb_trace_sigma[j] = qq_trace_factored_shared(
            K, BSigma_col_over_SV, K, KQU, have_KQU, pvj);
        pv_sigma[j] = pvj; }
    }

    _ss.stop();
    ngme_timing::Scope _sz(ngme_timing::rb_sec_Z_us());
    // Add Z-related RB trace: T = (dZ_j)^T A_i^T D A_i Z_i
    // where D is measurement precision (depends on noise settings)
    // Matches: T = dZ^T A^T D A Z
    // where D is measurement precision (depends on noise settings)
    // This whole section contributes exactly zero unless Z depends on theta.
    // Most operators build a Z that does not and dZ is then a full-size matrix
    // of ZEROS rather than an empty one, so the rows()/cols() guard below never
    // fired: every parameter paid for dZ^T (A^T D A) Z and a trace of the
    // result, to add zero. The fractional matern DOES carry a theta-dependent
    // Z (K and Z are the two halves of the rational approximation), so this has
    // to be decided from the values, not from the model.
    // NGME_RB_KEEP_ZERO_DZ=1 restores the old behaviour.
    static const bool keep_zero_dZ = [] {
      const char *e = std::getenv("NGME_RB_KEEP_ZERO_DZ");
      return e && *e && std::string(e) != "0";
    }();
    bool any_dZ = keep_zero_dZ;
    for (int j = 0; j < latents[i]->get_n_theta_K() && !any_dZ; ++j) {
      const auto &dZ_j = latents[i]->get_dZ(j);
      any_dZ = dZ_j.rows() > 0 && dZ_j.cols() > 0 && dZ_j.nonZeros() > 0;
    }

    // Precompute A_i^T D A_i once per latent -- and only if it will be used.
    SparseMatrix<double> ADA_i;
    if (any_dZ) {
      const auto &Ai = latents[i]->getA();
      if (!corr_measure) {
        // D = diag(1 / (sigma^2 V))
        VectorXd inv_noise_SV =
            noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
        ADA_i = Ai.transpose() * inv_noise_SV.asDiagonal() * Ai;
      } else {
        // Correlated case: D = Q_eps
        ADA_i = Ai.transpose() * Q_eps * Ai;
      }
    }

    if (!any_dZ) {
      // Skipping changes no VALUE but qq_trace() draws one probe seed per call on the
      // Hutchinson path, so skipping it silently would shift the random stream
      // and move every later draw. Advance the stream by exactly what the
      // skipped calls would have taken, so this optimisation is invisible in
      // the results rather than merely equivalent in expectation. (The exact
      // selected-inverse path draws nothing, hence the state test;
      // selinv_state_ is already decided here by the theta_K loop above.)
      if (selinv_state_ != 1)
        for (int j = 0; j < latents[i]->get_n_theta_K(); ++j)
          rng();
    }

    // Compute and accumulate per-parameter Z traces
    for (int j = 0; any_dZ && j < latents[i]->get_n_theta_K(); ++j) {
      const auto &dZ_j = latents[i]->get_dZ(j);
      if (dZ_j.rows() == 0 || dZ_j.cols() == 0 ||
          (!keep_zero_dZ && dZ_j.nonZeros() == 0))
        continue; // Z does not move with this parameter: the trace is zero
      const auto &Zi = latents[i]->getZ();
      // Local block in latent i's W-space
      SparseMatrix<double> Tloc = dZ_j.transpose() * ADA_i * Zi;
      // Embed into full W-space
      SparseMatrix<double> T(W_sizes, W_sizes);
      setSparseBlock(&T, woff, woff, Tloc);
      // Accumulate into rb_trace_K (same sign as K-term; Block compensates
      // later)
      { double pvj = 0.0; rb_trace_K[j] -= qq_trace(T, pvj); pv_K[j] += pvj; }
    }

    // Save per-latent RB traces at Block level for later gradient compensation
    rb_trace_K_latent[i] = rb_trace_K;
    rb_trace_sigma_latent[i] = rb_trace_sigma;
    rb_probe_var_K_latent[i] = pv_K;
    rb_probe_var_sigma_latent[i] = pv_sigma;
    n += latents[i]->get_V_size();
    woff += latents[i]->get_W_size();
  }

  // compute for theta_sigma
  ngme_timing::Scope _sn(ngme_timing::rb_sec_noise_us());
  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());
  // AZ does not depend on j, and get_AZ() already holds exactly this matrix,
  // rebuilt only when set_parameter_and_update() invalidates it. Assembling it
  // per j meant recomputing A_i Z_i for every theta_sigma on every Gibbs draw.
  const SparseMatrix<double> &AZ = get_AZ();
  const std::vector<int> noise_sig_free =
      free_sigma_cols(fix_theta_sigma_vec, (int)B_sigma.cols());
  for (int j = 0; j < n_theta_sigma; j++) {
    // The trace is added to the gradient of the j-th FREE theta_sigma, so it
    // has to be taken against that component's column of B_sigma.
    const int col_j = (j < (int)noise_sig_free.size()) ? noise_sig_free[j] : j;
    // T = (A Z)^T diag(B_sigma_j / noise_SV) (A Z), but never formed: AZ has
    // n_obs rows, so that triple product is by far the most expensive thing in
    // this function. trace_factored() applies the same three factors to the probe
    // block instead, where nothing bigger than n x N_iter is ever built.
    { double pvj = 0.0;
      rb_trace_noise_sigma[j] = qq_trace_factored(
          AZ, B_sigma.col(col_j).cwiseQuotient(noise_SV), AZ, pvj);
      if (rb_probe_var_noise_sigma.size() != n_theta_sigma)
        rb_probe_var_noise_sigma = VectorXd::Zero(n_theta_sigma);
      rb_probe_var_noise_sigma[j] = pvj; }
  }
  if (debug) {
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t_start)
                  .count();
    ngme_io::out() << "after compute trace (" << ms << " ms)" << std::endl;
  }
}

void BlockModel::assemble_dK() {
  int nrow = 0;
  int ncol = 0;
  for (int i = 0; i < n_latent; i++) {
    for (int j = 0; j < latents[i]->get_n_theta_K(); j++) {
      setSparseBlock(&block_dK[i][j], nrow, ncol, latents[i]->get_dK(j));
    }
    nrow += latents[i]->get_V_size();
    ncol += latents[i]->get_W_size();
  }
}

// generate output to R
Rcpp::List BlockModel::output() const {
  Rcpp::List latents_output;
  for (int i = 0; i < n_latent; i++) {
    latents_output.push_back((*latents[i]).output());
  }

  Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("noise") = Rcpp::List::create(
          Rcpp::Named("noise_type") = family,
          Rcpp::Named("theta_mu") = theta_mu,
          Rcpp::Named("theta_sigma") = theta_sigma,
          Rcpp::Named("theta_nu") = theta_nu, Rcpp::Named("V") = noise_V,
          Rcpp::Named("rho") = rho),
      Rcpp::Named("feff") = beta,
      // > 0 when the optimiser (hence the stored trajectory) used the
      // standardised NIG coordinates for the measurement noise.
      Rcpp::Named("merr_nig_std") = merr_nig_mode(),
      // Rcpp::Named("sampling_time")   = sampling_time.count(),
      Rcpp::Named("models") = latents_output
      // Rcpp::Named("log_likelihood")  = all_gaussian ? -log_likelihood() : 0
  );

  return out;
}

// AZ = [A_1 Z_1, ..., A_L Z_L]; depends on the parameters only, so it is
// cached and rebuilt by invalidate_AZ() when set_parameter_and_update() runs.
const SparseMatrix<double> &BlockModel::get_AZ() const {
  if (!AZ_valid || ngme_counters::cache_disabled()) {
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li = 0; li < n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    AZ_cached = AZ;
    AZ_valid = true;
  }
  return AZ_cached;
}

// H = sqrt(D) A Z, with D the measurement precision (and an extra sqrt_Rinv
// factor in the correlated case). Depends on AZ and on the measurement
// precision only, so it survives the latent V samplers untouched.
const SparseMatrix<double> &BlockModel::get_sqrt_AtSVA() const {
  if (sqrt_AtSVA_valid && !ngme_counters::cache_disabled())
    return sqrt_AtSVA_cached;
  const VectorXd sqrt_inv_noise_SV = meas_prec().cwiseSqrt();
  const SparseMatrix<double> &AZ = get_AZ();
  if (!corr_measure) {
    sqrt_AtSVA_cached = sqrt_inv_noise_SV.asDiagonal() * AZ;
  } else {
    sqrt_AtSVA_cached = sqrt_Rinv * sqrt_inv_noise_SV.asDiagonal() * AZ;
  }
  sqrt_AtSVA_valid = true;
  return sqrt_AtSVA_cached;
}

// G = sqrt(1/SV) K. Note this uses the *un-clamped* 1/SV that sampleW_VY()
// forms, not the floored one update_QQ() uses, so the caller passes it in.
const SparseMatrix<double> &BlockModel::get_G(const VectorXd &inv_SV) const {
  if (G_valid && !ngme_counters::cache_disabled())
    return G_cached;
  // As in get_sqrt_AtSVA(), evaluate the square root eagerly rather than
  // letting asDiagonal() hold a lazy expression across the sparse product.
  const VectorXd sqrt_inv_SV = inv_SV.cwiseSqrt();
  G_cached = sqrt_inv_SV.asDiagonal() * K;
  G_valid = true;
  return G_cached;
}

// Measurement block of QQ: H'H in the uncorrelated case (H as above), and
// (AZ)' Q_eps (AZ) in the correlated one. Same dependencies as H.
const SparseMatrix<double> &BlockModel::get_QQ_measure() const {
  if (QQ_measure_valid && !ngme_counters::cache_disabled())
    return QQ_measure;
  ngme_timing::Scope _s(ngme_timing::qq_measure_us());
  if (!corr_measure) {
    const SparseMatrix<double> &H = get_sqrt_AtSVA();
    QQ_measure = H.transpose() * H;
  } else {
    const SparseMatrix<double> &AZ = get_AZ();
    QQ_measure = AZ.transpose() * Q_eps * AZ;
  }
  QQ_measure_valid = true;
  return QQ_measure;
}

// update Q_eps, dQ_eps, and compute trace as sum_ij dQ_ij * Q^-1_ij
void BlockModel::update_Q_eps(double rho) {
  invalidate_measurement();
  // update Q_eps
  for (int i = 0; i < Q_eps.outerSize(); i++) {
    for (SparseMatrix<double>::InnerIterator it(Q_eps, i); it; ++it) {
      if (it.row() == it.col()) {
        int idx = it.row();
        it.valueRef() = 1.0 / (pow(noise_sigma(idx), 2) * noise_V(idx));
        if (has_correlation[idx])
          it.valueRef() /= (1 - rho * rho);
      } else {
        double tmp = noise_sigma(it.row()) * noise_sigma(it.col()) *
                     sqrt(noise_V(it.row()) * noise_V(it.col()));
        it.valueRef() = -rho / ((1 - rho * rho) * tmp);
      }
    }
  }

  // update sqrt_Rinv
  for (int i = 0; i < sqrt_Rinv.outerSize(); i++) {
    for (SparseMatrix<double>::InnerIterator it(sqrt_Rinv, i); it; ++it) {
      if (it.row() == it.col()) {
        int idx = it.row();
        it.valueRef() = has_correlation[idx] ? (sqrt(1 + rho) + sqrt(1 - rho)) /
                                                   (2 * sqrt(1 - rho * rho))
                                             : 1;
      } else {
        it.valueRef() =
            (-sqrt(1 + rho) + sqrt(1 - rho)) / (2 * sqrt(1 - rho * rho));
      }
    }
  }

  // update dQ_eps, and compute trace as sum_ij dQ_ij * Q^-1_ij
  for (int i = 0; i < dQ_eps.outerSize(); i++) {
    for (SparseMatrix<double>::InnerIterator it(dQ_eps, i); it; ++it) {
      if (it.row() == it.col()) {
        int idx = it.row();
        double tmp = pow((1 - rho * rho) * noise_sigma(idx), 2) * noise_V(idx);
        it.valueRef() = 2.0 * rho / tmp;
      } else {
        int r = it.row();
        int c = it.col();
        double tmp = pow(1 - rho * rho, 2) * noise_sigma(r) * noise_sigma(c) *
                     sqrt(noise_V(r) * noise_V(c));
        it.valueRef() = -(1 + rho * rho) / tmp;
      }
    }
  }

  // compute logdet of Q_eps (not needed)
  // lhs = Q_eps.logdet();
  // VectorXd SV =
  // noise_sigma.array().pow(2).matrix().cwiseProduct(noise_prevV); double lhs =
  // 0; int i = 0; while (i < n_obs) {
  //   if (has_correlation[i]) {
  //     lhs += -log((1 - rho * rho) * SV[i] * SV[i + 1]);
  //     i += 2;
  //   } else {
  //     lhs += -log(SV[i]);
  //     i += 1;
  //   }
  // }
  // VectorXd res = get_residual();
  // double rhs = res.dot(Q_eps * res);
}

// Structural signature of the QQ that chol_QQ.analyze() last saw. Comparing it
// is O(nnz) and far cheaper than redoing the fill-reducing ordering and the
// symbolic factorization.
void BlockModel::record_QQ_pattern() {
  const int outer = QQ.outerSize() + 1;
  const int nnz = static_cast<int>(QQ.nonZeros());
  QQ_pat_outer.assign(QQ.outerIndexPtr(), QQ.outerIndexPtr() + outer);
  QQ_pat_inner.assign(QQ.innerIndexPtr(), QQ.innerIndexPtr() + nnz);
  QQ_analyzed = true;
}

bool BlockModel::QQ_pattern_changed() const {
  if (!QQ_analyzed)
    return true;
  if (QQ_pat_outer.size() != static_cast<size_t>(QQ.outerSize()) + 1)
    return true;
  if (QQ_pat_inner.size() != static_cast<size_t>(QQ.nonZeros()))
    return true;
  if (!std::equal(QQ_pat_outer.begin(), QQ_pat_outer.end(),
                  QQ.outerIndexPtr()))
    return true;
  return !std::equal(QQ_pat_inner.begin(), QQ_pat_inner.end(),
                     QQ.innerIndexPtr());
}

namespace {
// Diagonal nudge used to rescue a QQ Cholesky that failed on a transient SGD
// excursion; relative to the largest diagonal entry of QQ, escalating by a
// factor of ten per attempt (1e-10 up to 1e-5 of the diagonal scale).
constexpr double QQ_JITTER_BASE = 1e-10;
constexpr int QQ_JITTER_ATTEMPTS = 6;
} // namespace

// Check that QQ's pattern really is the union of Q's and the measurement
// block's, and record the nnz the check was made against. Both operands are
// CSC with sorted inner indices, so one merge walk per column visits every
// entry; anything left unplaced means the fast merge in update_QQ() would
// silently drop terms, and the guard nnz are reset so it is not taken.
void BlockModel::record_qq_add_pattern(const SparseMatrix<double> &Qm,
                                       const SparseMatrix<double> &Me) {
  const int *Qp = Qm.outerIndexPtr(), *Qi = Qm.innerIndexPtr();
  const int *Mp = Me.outerIndexPtr(), *Mi = Me.innerIndexPtr();
  const int *Tp = QQ.outerIndexPtr(), *Ti = QQ.innerIndexPtr();
  for (int c = 0; c < QQ.outerSize(); ++c) {
    int qp = Qp[c], mp = Mp[c];
    const int qe = Qp[c + 1], me = Mp[c + 1];
    for (int t = Tp[c]; t < Tp[c + 1]; ++t) {
      const int row = Ti[t];
      if (qp < qe && Qi[qp] == row) ++qp;
      if (mp < me && Mi[mp] == row) ++mp;
    }
    if (qp != qe || mp != me) {
      qq_map_q_nnz_ = qq_map_meas_nnz_ = -1;
      return;
    }
  }
  qq_map_q_nnz_ = (long long)Qm.nonZeros();
  qq_map_meas_nnz_ = (long long)Me.nonZeros();
}

void BlockModel::update_QQ() {
  ngme_counters::bump(ngme_counters::QQ_builds);
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
  inv_SV = inv_SV.cwiseMax(1e-8);

  // update Q and QQ. Only the latent block K' diag(1/SV) K has to be redone on
  // every draw; the measurement block is cached and reused.
  { ngme_timing::Scope _s(ngme_timing::qq_assemble_us());
    { ngme_timing::Scope _p(ngme_timing::qq_prod_us());
      // Caching K^T across Gibbs draws was tried here (K is constant over a
      // sweep while 1/SV is not) and measured no better. What DOES pay is
      // caching the destination slots of the triple product itself: K's
      // pattern is fixed for the fit, so Eigen's symbolic phase is repeated
      // every draw for nothing. The cache declines and this falls back to the
      // direct product when the scatter list would be too large.
      if (!qq_ata_.refill(K, inv_SV, Q, qq_ata_budget_))
        Q = K.transpose() * inv_SV.asDiagonal() * K; }
    { ngme_timing::Scope _a(ngme_timing::qq_add_us());
      const SparseMatrix<double> &Me = get_QQ_measure();
      // QQ's pattern is the union of the two operands' and does not move while
      // they keep their nnz, so it is established once and thereafter only the
      // values are refilled -- which matters because this runs once per Gibbs
      // draw of a non-Gaussian model.
      const bool pattern_ok = qq_map_q_nnz_ == (long long)Q.nonZeros() &&
                              qq_map_meas_nnz_ == (long long)Me.nonZeros() &&
                              QQ.rows() == Q.rows() && QQ.isCompressed() &&
                              !ngme_counters::cache_disabled();
      if (!pattern_ok) {
        QQ = Q + Me;
        QQ.makeCompressed();
        record_qq_add_pattern(Q, Me);
      } else {
        // Merge both operands straight into QQ's existing storage. All three
        // are CSC with sorted inner indices, so one walk per column fills it
        // with no allocation.
        const int *Qp = Q.outerIndexPtr(), *Qi = Q.innerIndexPtr();
        const int *Mp = Me.outerIndexPtr(), *Mi = Me.innerIndexPtr();
        const int *Tp = QQ.outerIndexPtr(), *Ti = QQ.innerIndexPtr();
        const double *qv = Q.valuePtr(), *mv = Me.valuePtr();
        double *v = QQ.valuePtr();
        for (int c = 0; c < QQ.outerSize(); ++c) {
          int qp = Qp[c], mp = Mp[c];
          const int qe = Qp[c + 1], me = Mp[c + 1];
          for (int t = Tp[c]; t < Tp[c + 1]; ++t) {
            const int row = Ti[t];
            double acc = 0.0;
            if (qp < qe && Qi[qp] == row) acc += qv[qp++];
            if (mp < me && Mi[mp] == row) acc += mv[mp++];
            v[t] = acc;
          }
        }
      } } }
  if (robust) {
    QQ = 0.5 * (QQ + SparseMatrix<double>(QQ.transpose()));
    // double mean_diag = QQ.diagonal().mean();
    // double jitter = std::max(1e-8, 1e-3 * mean_diag);
    double jitter = 1e-8;
    QQ.diagonal().array() += jitter;
  }
  // Pack, re-run the symbolic phase if the pattern moved, and factorize.
  // Storage must be packed before the pattern can be compared (and before the
  // solvers see it); that is a no-op unless the robust branch above had to
  // insert a missing diagonal entry. The symbolic phase is only needed when
  // the sparsity pattern really moved. It does for rational approximations only.
  // For every other model the pattern is fixed and analyze() runs exactly once.
  auto factorize = [this]() {
    if (!QQ.isCompressed())
      QQ.makeCompressed();
    if (QQ_pattern_changed()) {
      ngme_counters::bump(ngme_counters::QQ_analyzes);
      { ngme_timing::Scope _s(ngme_timing::qq_symbolic_us()); chol_QQ.analyze(QQ); }
      setup_qq_probing();
      record_QQ_pattern();
    }
    { ngme_timing::Scope _s(ngme_timing::qq_numeric_us()); chol_QQ.compute(QQ); }
    return chol_QQ.factorization_success();
  };

  QQ_valid = factorize();
  if (QQ_valid)
    return;

  // A Cholesky failure here is almost always a transient SGD excursion rather
  // than a genuinely indefinite QQ: the diagonal is still strictly positive
  // and the asymmetry sits at round-off level, but the iterate has drifted
  // close enough to singular that the factorization gives up. Nudging the
  // diagonal recovers it, and a perturbation this small is orders of magnitude
  // below the Monte Carlo error of the sampled gradient, so escalate the nudge
  // a few times before declaring the fit dead. Aborting instead would throw
  // away a run that the next iterate usually recovers from on its own.
  const double scale = std::max(1.0, QQ.diagonal().cwiseAbs().maxCoeff());
  SparseMatrix<double> Id(QQ.rows(), QQ.cols());
  Id.setIdentity();
  double applied = 0.0;
  for (int attempt = 0; attempt < QQ_JITTER_ATTEMPTS && !QQ_valid; ++attempt) {
    // QQ_JITTER_BASE, 10x that, 100x that, ... relative to the largest
    // diagonal entry. Going through the identity rather than
    // QQ.diagonal().array() += ... inserts any structurally absent diagonal
    // entry instead of silently skipping it.
    const double target = scale * QQ_JITTER_BASE * pow(10.0, attempt);
    SparseMatrix<double> jittered = QQ + (target - applied) * Id;
    QQ.swap(jittered);
    applied = target;
    QQ_valid = factorize();
  }
  if (QQ_valid)
    return;

  // Out of options. Report the diagnostics through the exception instead of
  // printing them here: update_QQ() runs inside the OpenMP region of the
  // parallel SGD loop, and Rcpp's output streams are R API calls that must not
  // be made off the main thread. estimate.cpp captures the message per chain
  // and re-raises it on the main thread.
  std::ostringstream msg;
  msg << "Measurement precision QQ is not SPD at iteration " << curr_iter
      << " (smallest diagonal entry " << QQ.diagonal().minCoeff()
      << ", asymmetry " << (QQ - SparseMatrix<double>(QQ.transpose())).norm()
      << "); adding up to " << applied
      << " to the diagonal did not restore positive definiteness. The model is "
         "probably unidentifiable at these parameter values, or the optimizer "
         "step size is too large.";
  throw std::runtime_error(msg.str());
}
