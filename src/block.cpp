// Implementation for block model and block_rep

#include "block.h"
#include "prior.h"
#include "sample_rGIG.h"
#include <random>
#include <cmath>
#include <iterator>
#include <chrono>

using std::pow;

// -------------- Block Model class ----------------
BlockModel::BlockModel(
  const Rcpp::List& block_model,
  unsigned long seed
) :
  rng               (seed),
  X                 (Rcpp::as<MatrixXd>      (block_model["X"])),
  Y                 (Rcpp::as<VectorXd>      (block_model["Y"])),
  W_sizes           (Rcpp::as<int>           (block_model["W_sizes"])),
  V_sizes           (Rcpp::as<int>           (block_model["V_sizes"])),
  beta              (Rcpp::as<VectorXd>      (block_model["feff"])),
  n_obs             (Y.size()),
  n_la_params       (Rcpp::as<int>           (block_model["n_la_params"])),
  n_feff            (beta.size()),
  n_merr            (Rcpp::as<int>           (block_model["n_merr"])),
  n_repl            (Rcpp::as<int>           (block_model["n_repl"])),
  // corr_measure      (Rcpp::as<bool>          (block_model["corr_measure"])),
  // cor_rows          (),
  // cor_cols          (),
  // has_correlation   (),
  // n_corr_pairs      (0),
  Q_eps             (n_obs, n_obs),
  dQ_eps            (n_obs, n_obs),
  n_params          (n_la_params + n_feff + n_merr),

  debug             (false),
  A                 (n_obs, W_sizes),
  K                 (V_sizes, W_sizes),
  Q                 (W_sizes, W_sizes),
  QQ                (W_sizes, W_sizes),

  // var               (Var(Rcpp::as<Rcpp::List> (block_model["noise"]), rng())),
  p_vec             (n_obs),
  a_vec             (n_obs),
  b_vec             (n_obs),
  noise_V           (VectorXd::Ones(n_obs)),
  noise_prevV       (VectorXd::Ones(n_obs)),

  curr_iter         (0),
  all_gaussian     (Rcpp::as<bool>       (block_model["all_gaussian"])),
  par_string       (Rcpp::as<string>     (block_model["par_string"]))
{
  // 1. Init controls
  Rcpp::List control_ngme = block_model["control_ngme"];
    // const double stepsize = control_ngme["stepsize"];
    bool init_sample_W = Rcpp::as<bool> (control_ngme["init_sample_W"]);
    n_gibbs     =  Rcpp::as<int>    (control_ngme["n_gibbs_samples"]);
    debug       = Rcpp::as<bool>   (control_ngme["debug"]);
    rao_blackwell = Rcpp::as<bool> (control_ngme["rao_blackwellization"]);
    int n_trace_iter = Rcpp::as<int> (control_ngme["n_trace_iter"]);
    int solver_type = control_ngme.containsElementNamed("solver_type") ? Rcpp::as<int> (control_ngme["solver_type"]) : 0;

    // reduce_var    =  Rcpp::as<bool>   (control_ngme["reduce_var"]);
    // reduce_power  =  Rcpp::as<double> (control_ngme["reduce_power"]);
    // threshold   =  Rcpp::as<double> (control_ngme["threshold"]);

if (debug) std::cout << "Begin Block Constructor" << std::endl;

  // 2. Init Fixed effects
  fix_flag[block_fix_beta]   = Rcpp::as<bool>        (control_ngme["fix_feff"]);
  if (beta.size() == 0) fix_flag[block_fix_beta]  = true;

  // 4. Init latent models
  Rcpp::List latents_in = block_model["models"];
  n_latent = latents_in.size(); // how many latent model
  if (n_latent == 0) rao_blackwell = false;
  for (int i=0; i < n_latent; ++i) {
    // construct acoording to models
    Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);
    latent_in["solver_type"] = solver_type;
    latent_in["n_trace_iter"] = n_trace_iter;
    unsigned long latent_seed = seed + (i+1) * 1000;
    latents.push_back(std::make_shared<Latent>(latent_in, latent_seed));
  }

if (debug) std::cout << "before set block A" << std::endl;
  /* Init A */
  int n = 0;
  for (std::vector<std::shared_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
    setSparseBlock(&A, 0, n, (*it)->getA());
    n += (*it)->get_W_size();
  }
if (debug) std::cout << "After set block K" << std::endl;

  // 5. Init measurement noise (consider corr_measure)
  Rcpp::List noise_in   = block_model["noise"];
    fix_flag[block_fix_theta_mu]     = Rcpp::as<bool> (noise_in["fix_theta_mu"]);
    
    // Handle vector-based fix_theta_sigma - simplified for block model
    Rcpp::LogicalVector fix_theta_sigma_r = Rcpp::as<Rcpp::LogicalVector>(noise_in["fix_theta_sigma"]);
    fix_flag[block_fix_theta_sigma] = std::all_of(fix_theta_sigma_r.begin(), fix_theta_sigma_r.end(), [](bool x) { return x; });

    fix_flag[blcok_fix_V]            = Rcpp::as<bool> (noise_in["fix_V"]);
    fix_flag[block_fix_theta_nu]     = noise_in.containsElementNamed("fix_theta_nu") ? Rcpp::as<bool> (noise_in["fix_theta_nu"]) : false;
    fix_flag[block_fix_rho]          = noise_in.containsElementNamed("fix_rho") ? Rcpp::as<bool> (noise_in["fix_rho"]) : false;
    shared_sigma = noise_in.containsElementNamed("shared_sigma") ? Rcpp::as<bool> (noise_in["shared_sigma"]) : false;

    B_mu          = (Rcpp::as<MatrixXd>      (noise_in["B_mu"])),
    theta_mu      = (Rcpp::as<VectorXd>      (noise_in["theta_mu"])),
    n_theta_mu    = (Rcpp::as<int>           (noise_in["n_theta_mu"])),

    B_sigma       = (Rcpp::as<MatrixXd>      (noise_in["B_sigma"])),
    theta_sigma   = (Rcpp::as<VectorXd>      (noise_in["theta_sigma"])),
    n_theta_sigma = (Rcpp::as<int>           (noise_in["n_theta_sigma"])),

    B_nu          = (Rcpp::as<MatrixXd>      (noise_in["B_nu"])),
    theta_nu      = (Rcpp::as<VectorXd>      (noise_in["theta_nu"])),
    n_theta_nu    = (Rcpp::as<int>           (noise_in["n_theta_nu"])),

    rb_trace_noise_sigma = VectorXd::Zero(n_theta_sigma),

    family = Rcpp::as<string>  (noise_in["noise_type"]);
    noise_mu    = B_mu * theta_mu;
    noise_sigma = (B_sigma * theta_sigma).array().exp();
    noise_nu    = (B_nu * theta_nu).array().exp();

    rho   = Rcpp::as<VectorXd> (noise_in["rho"]); 
    n_rho = noise_in.containsElementNamed("n_rho") ? Rcpp::as<int>      (noise_in["n_rho"]) : 0;
    nu_lower_bound = noise_in.containsElementNamed("nu_lower_bound") ? Rcpp::as<double> (noise_in["nu_lower_bound"]) : 0.0;
    corr_measure = Rcpp::as<bool> (noise_in["corr_measurement"]);

    // init priors for noise_parameter
    Rcpp::List prior_list = Rcpp::as<Rcpp::List> (noise_in["prior_mu"]);
        prior_mu_type  = Rcpp::as<string> (prior_list[0]);
        prior_mu_param = Rcpp::as<VectorXd> (prior_list["param"]);
    prior_list = Rcpp::as<Rcpp::List> (noise_in["prior_sigma"]);
        prior_sigma_type  = Rcpp::as<string> (prior_list["type"]);
        prior_sigma_param = Rcpp::as<VectorXd> (prior_list["param"]);
    prior_list = Rcpp::as<Rcpp::List> (noise_in["prior_nu"]);
        prior_nu_type  = Rcpp::as<string> (prior_list["type"]);
        prior_nu_param = Rcpp::as<VectorXd> (prior_list["param"]);

    if (family != "normal") {
      NoiseUtil::update_gig(family, noise_nu, p_vec, a_vec, b_vec);
    }

    if (corr_measure) {
      cor_cols = Rcpp::as<vector<int>> (noise_in["cor_cols"]);
      cor_rows = Rcpp::as<vector<int>> (noise_in["cor_rows"]);
      has_correlation = Rcpp::as<vector<bool>> (noise_in["has_correlation"]);

// print has_correlation
// std::cout << "has_correlation = ";
// for (int i=0; i < has_correlation.size(); i++) {
//   std::cout << has_correlation[i] << " ";
// }

      n_corr_pairs = Rcpp::as<int> (noise_in["n_corr_pairs"]);
      vector<Triplet<double>> Q_eps_triplet, dQ_eps_triplet;
      for (int i=0; i < cor_cols.size(); ++i) {
        Q_eps_triplet.push_back(Triplet<double>(cor_rows[i], cor_cols[i], cor_rows[i] == cor_cols[i]));
        if (has_correlation[cor_rows[i]]) {
          // ignore uncorrelated locations
          dQ_eps_triplet.push_back(Triplet<double>(cor_rows[i], cor_cols[i], cor_rows[i] == cor_cols[i]));
        }
      }
      SparseMatrix<double> Q_eps_lower (n_obs, n_obs);
      Q_eps_lower.setFromTriplets(Q_eps_triplet.begin(), Q_eps_triplet.end());
      Q_eps = Q_eps_lower.selfadjointView<Lower>();
      sqrt_Rinv = Q_eps; // initialize sqrt_Rinv as Q_eps

      SparseMatrix<double> dQ_lower (n_obs, n_obs);
      dQ_lower.setFromTriplets(dQ_eps_triplet.begin(), dQ_eps_triplet.end());
      dQ_eps = dQ_lower.selfadjointView<Lower>();

      update_Q_eps(rho(0));
// std::cout << "rho(0) = " << rho(0) << std::endl;
// std::cout << "noise_sigma = " << noise_sigma.transpose() << std::endl;
// std::cout << "Init Q_eps: \n" << Q_eps << std::endl;
    }

if (debug) std::cout << "After block construct noise" << std::endl;

  // 6. Fix V and init V
  if (noise_in.containsElementNamed("V") && !Rf_isNull(noise_in["V"])) {
      noise_V = Rcpp::as< VectorXd > (noise_in["V"]);
      noise_prevV = noise_V;
  }

  // 7. Init solvers
  assemble();
if (debug) std::cout << "After assemble" << std::endl;

  if (n_latent > 0) {
    VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
    Q = K.transpose() * inv_SV.asDiagonal() * K;
    // Build AZ for measurement term
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col0 = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col0, AiZi);
      col0 += latents[li]->get_W_size();
    }
    if (!corr_measure && shared_sigma) {
      QQ = Q + AZ.transpose() * noise_V.cwiseInverse().asDiagonal() * AZ;
    } else if (!corr_measure) {
      QQ = Q + AZ.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * AZ;
    } else{
      QQ = Q + AZ.transpose() * Q_eps * AZ;
    }

    // Initialize solver with requested backend and Hutchinson iters; QQ is SPD
    chol_QQ.init(QQ.rows(), n_trace_iter, /*symmetric*/ true, solver_type);
    chol_QQ.analyze(QQ);
    chol_QQ.compute(QQ);
  }

  // 8. optimizer related
  // stepsizes = VectorXd::Constant(n_params, stepsize);
  steps_to_threshold = VectorXd::Constant(n_params, 0);
  indicate_threshold = VectorXd::Constant(n_params, 0);

// if (debug) std::cout << "After init solver && before sampleW_V" << std::endl;

  // Default: enable RB for all-Gaussian models and use conditional W
  if (rao_blackwell) {
    // Initialize block_dKs of length n_latent
    block_dK.resize(n_latent);
    for (int i=0; i < n_latent; i++) {
      block_dK[i].resize(latents[i]->get_n_theta_K());
      for (int j=0; j < latents[i]->get_n_theta_K(); j++) {
        block_dK[i][j] = SparseMatrix<double>(V_sizes, W_sizes);
      }
    }
  }
  
  // Initialize gradient covariance matrix
  grad_covariance = MatrixXd::Zero(n_params, n_params);
  
if (debug) std::cout << "End Block Constructor" << std::endl;
}

  void BlockModel::burn_in(int iterations) {
      for (int i=0; i < iterations; i++) {
          sample_cond_V(false);
          sampleW_VY(true);
          sample_cond_noise_V();
      }
  }

void BlockModel::setW(const VectorXd& W) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->setW(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::set_cond_W(const VectorXd& W) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->set_cond_W(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::setPrevW(const VectorXd& W) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->setPrevW(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::setPrevV(const VectorXd& V) {
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_V_size();
    (*it)->setPrevV(V.segment(pos, size));
    pos += size;
  }
}

// sample W|VY
void BlockModel::sampleW_VY(bool burn_in) {
  if (n_latent==0) return;
  // Ensure QQ is consistent with current K, Z, and measurement precision before sampling
  update_QQ();
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());

  // M = K' * inv(SV) * mean + Z'^ A'^ inv(Sigma) * (Y - X * beta - (1 - V) mu)
  VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean();

  if (!corr_measure && shared_sigma) {
    // special case of shared_sigma, replace Y by Y/sigma_e
    // add Z'^A' * (Y/sigma - X beta - (1-V)mu)
    VectorXd rhs = get_residual_part() - Y + Y.cwiseQuotient(noise_sigma);
    // rhs = Y/sigma - X beta - (1-V)mu
    // Build AZ and apply
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    M += AZ.transpose() * rhs;
  }
  else if (!corr_measure) {
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    M += AZ.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * get_residual_part();
  } else {
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    M += AZ.transpose() * Q_eps * get_residual_part();
  }

  SparseMatrix<double> G = inv_SV.cwiseSqrt().asDiagonal() * K;
  SparseMatrix<double> H = get_sqrt_AtSVA();
  unsigned long seed1 = rng();
  unsigned long seed2 = rng();
  VectorXd z1 = NoiseUtil::rnorm_vec(G.rows(), 0, 1, seed1);
  VectorXd z2 = NoiseUtil::rnorm_vec(H.rows(), 0, 1, seed2);

  // sample W ~ N(QQ^-1*M, QQ^-1)
  
  // Sampling method using cholesky decomposition, using matrixL()
  // W = QQ^-1*M + L^-1*z1 + L^-T*z2
  // VectorXd W = chol_QQ.rMVN(M, z1);

  // Sampling method using tricks, purely solve, not matrixL()
  VectorXd W_draw = chol_QQ.rMVN(G, H, M, z1, z2);
  setW(W_draw);

  if (rao_blackwell && !burn_in) {
    VectorXd W_mean = chol_QQ.solve(M);
    set_cond_W(W_mean);
  }
// std::cout << "size of W and time of sampling is " << W.size() << " " << time << std::endl;
// if (debug) std::cout << "Finish sampling W" << std::endl;
}

// ---------------- get, set update gradient ------------------
// order is Latent, merr, feff
VectorXd BlockModel::get_parameter() {
if (debug) std::cout << "Start block get_parameter"<< std::endl;
    VectorXd thetas (n_params);
    int pos = 0;
    for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      VectorXd theta = (*it)->get_parameter();
      thetas.segment(pos, theta.size()) = theta;
      pos += theta.size();
    }
    thetas.segment(n_la_params, n_merr) = get_theta_merr();

    if (!fix_flag[block_fix_beta]) 
      thetas.segment(n_la_params + n_merr, n_feff) = beta;
    
if (debug) std::cout << "Finish block get_parameter"<< std::endl;
    return thetas;
}

void BlockModel::set_parameter(const VectorXd& Theta) {
// if (debug) std::cout << "Start block set_parameter"<< std::endl;
// std::chrono::steady_clock::time_point startTime, endTime; startTime = std::chrono::steady_clock::now();
  int pos = 0;
  for (std::vector<std::shared_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
    int theta_len = (*it)->get_n_params();
    VectorXd theta = Theta.segment(pos, theta_len);
    (*it)->set_parameter(theta, rao_blackwell);
    pos += theta_len;
  }

  // measurement noise
  set_theta_merr(Theta.segment(n_la_params, n_merr));
  if (family!="normal") NoiseUtil::update_gig(family, noise_nu, p_vec, a_vec, b_vec);

  // fixed effects
  if (!fix_flag[block_fix_beta]) {
    beta = Theta.segment(n_la_params + n_merr, n_feff);
  }

  assemble(); //update K,dK,d2K after
// endTime = std::chrono::steady_clock::now(); update_time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime); std::cout << "block set time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << std::endl;
  curr_iter++;
  if (n_latent > 0 && family == "normal") { update_QQ(); }
}

// --------- Fiexed effects and Measurement Error ---------------
VectorXd BlockModel::grad_beta() {
  VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());

  VectorXd residual = get_residual(rao_blackwell); // + X * beta;
  VectorXd grads = X.transpose() * noise_inv_SV.asDiagonal() * residual;

// special case
if (shared_sigma) {
  residual = residual - Y + Y.cwiseQuotient(noise_sigma); // replace Y by Y/sigma_e
  grads = X.transpose() * residual;
}
  //  * residual.cwiseQuotient(noise_sigma);

  // MatrixXd hess = X.transpose() * noise_inv_SV.asDiagonal() * X;
  // grads = grads / A.rows();
  // grads = hess.ldlt().solve(grads);

  // std::cout << "grads of beta=" << grads << std::endl;
    return -grads;
}

VectorXd BlockModel::grad_theta_mu() {
  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());

  if (shared_sigma) {noise_SV = noise_V;}
  VectorXd residual = get_residual(rao_blackwell);
  VectorXd grad = VectorXd::Zero(n_theta_mu);
  for (int l=0; l < n_theta_mu; l++) {
      grad(l) = (noise_V - VectorXd::Ones(n_obs)).cwiseProduct(B_mu.col(l).cwiseQuotient(noise_SV)).dot(residual);
      // add prior
      // grad(l) += PriorUtil::d_log_dens(prior_mu_type, prior_mu_param, theta_mu(l));
  }
  return -grad;
}

VectorXd BlockModel::grad_theta_sigma() {
  VectorXd grad = VectorXd::Zero(n_theta_sigma);
  VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);

  VectorXd residual = get_residual(rao_blackwell);
  VectorXd vsq = (residual).array().pow(2).matrix().cwiseQuotient(noise_SV);
  VectorXd tmp1 = vsq - VectorXd::Ones(n_obs);
  grad = B_sigma.transpose() * tmp1; // dℓ/dθ for generic case

  // Special case: Y|W ~ N(sigma_e (Xβ + A Z W), sigma_e)
  // dsigma_e/d theta contribution in log-lik form under log-σ paramization
  if (shared_sigma) {
    VectorXd mean = - get_residual(false) + Y;  // AW+Xβ
    VectorXd Y_by_sigma = Y.cwiseQuotient(noise_sigma);
    VectorXd Y_by_sigma_sq = Y_by_sigma.array().pow(2);
    grad = B_sigma.transpose() * (- VectorXd::Ones(n_obs) + Y_by_sigma_sq - Y_by_sigma.cwiseProduct(mean));
  }

  // Rao-Blackwell term from log|QQ|: keep sign consistent with theta_K path
  // theta_K: grad_accum -= rb_trace_K; then return -grad_accum
  // Here: do grad -= rb_trace_noise_sigma; then return -grad
  if (rao_blackwell) grad += rb_trace_noise_sigma;

  // add prior (disabled by default)
  // for (int l=0; l < n_theta_sigma; l++) {
  //   grad(l) += PriorUtil::d_log_dens(prior_sigma_type, prior_sigma_param, theta_sigma(l));
  // }

  // Return gradient of negative log-likelihood
  return -grad;
}

VectorXd BlockModel::get_theta_merr() const {
  VectorXd theta_merr = VectorXd::Zero(n_merr);

  if (!fix_flag[block_fix_theta_mu])
    theta_merr.segment(0, n_theta_mu) = theta_mu;
  if (!fix_flag[block_fix_theta_sigma])
    theta_merr.segment(n_theta_mu, n_theta_sigma) = theta_sigma;
  if (!fix_flag[block_fix_theta_nu])
    theta_merr.segment(n_theta_mu + n_theta_sigma, n_theta_nu) = theta_nu;
  
  if (corr_measure && !fix_flag[block_fix_rho])
    theta_merr(n_merr-1) = rho2th(rho(0));

  return theta_merr;
}

VectorXd BlockModel::grad_theta_merr() {
  VectorXd grad = VectorXd::Zero(n_merr);

  if (!fix_flag[block_fix_theta_mu])     grad.segment(0, n_theta_mu) = grad_theta_mu();
  if (!fix_flag[block_fix_theta_sigma])  grad.segment(n_theta_mu, n_theta_sigma) = grad_theta_sigma();
  if (!fix_flag[block_fix_theta_nu]) {
    grad.segment(n_theta_mu + n_theta_sigma, n_theta_nu) = NoiseUtil::grad_theta_nu(family, B_nu, noise_nu, noise_V, noise_prevV);
    // add prior
    // grad(n_theta_mu + n_theta_sigma) -= PriorUtil::d_log_dens(prior_nu_type, prior_nu_param, noise_nu);
  }

  // grad of theta_rho
  if (corr_measure && !fix_flag[block_fix_rho]) {
    // Q_eps_solver.factorize(Q_eps);
    double trace = 0.5 * 2 * rho(0)/(1-rho(0)*rho(0)) * n_corr_pairs;
    VectorXd res = get_residual();
    double drhs = -0.5 * (res).dot(dQ_eps * res);
    grad(n_merr-1) = trace + drhs;
    grad(n_merr-1) *= -dtheta_th(rho(0));
// std::cout << "drhs = " << drhs << std::endl;
// std::cout << "trace = " << trace << std::endl;
// std::cout << "grad of rho=" << grad(n_merr-1) << std::endl;
  }

  return grad;
}

void BlockModel::set_theta_merr(const VectorXd& theta_merr) {
// if (debug) std::cout << "start set theta_merr" << std::endl;
  if (!fix_flag[block_fix_theta_mu])     
    theta_mu = theta_merr.segment(0, n_theta_mu);
  if (!fix_flag[block_fix_theta_sigma])  
    theta_sigma = theta_merr.segment(n_theta_mu, n_theta_sigma);
  if (!fix_flag[block_fix_theta_nu]) 
    theta_nu = theta_merr.segment(n_theta_mu + n_theta_sigma, n_theta_nu);

  if (family != "normal" && theta_nu(0) > log(1e4)) theta_nu(0) = 1e4;

  // update rho, and Q_eps
  if (corr_measure && !fix_flag[block_fix_rho]) {
    rho(0) = th2rho(theta_merr(n_merr-1));
    update_Q_eps(rho(0));
  }

  // update mu, sigma
  noise_mu    = (B_mu * theta_mu);
  noise_sigma = (B_sigma * theta_sigma).array().exp();
  noise_nu    = (B_nu * theta_nu).array().exp();
  noise_nu    = noise_nu.cwiseMax(nu_lower_bound);

// show the Q construction
// std::cout << "Q_eps == \n" << Q_eps << std::endl;
// std::cout << "dQ_eps == \n" << dQ_eps << std::endl;
}

void BlockModel::sample_cond_noise_V(bool posterior) {
  if (family == "normal" || fix_flag[blcok_fix_V]) return;
  noise_prevV = noise_V;

  if (posterior) {
    VectorXd a_inc_vec = noise_mu.cwiseQuotient(noise_sigma).array().pow(2);
    VectorXd b_inc_vec = (get_residual() + noise_V.cwiseProduct(noise_mu)).cwiseQuotient(noise_sigma).array().pow(2);
    // VectorXd b_inc_vec = (Y - A * getW() - X * beta - noise_mu).cwiseQuotient(noise_sigma).array().pow(2);
    VectorXd a_vec_new = a_vec + a_inc_vec;
    VectorXd b_vec_new = b_vec + b_inc_vec;

    if (!corr_measure) {
      double dim = 1;
      VectorXd p_vec_new = p_vec - VectorXd::Constant(n_obs, 0.5 * dim);
      // noise_V = rGIG_cpp(p_vec_new, a_vec_new, b_vec_new, rng());
      NoiseUtil::sample_V(noise_V, family, p_vec_new, a_vec_new, b_vec_new, rng);
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
          noise_V[i] = noise_V[i+1] = rGIG_cpp(p_vec[i], a_vec[i], b_vec[i], rng());
          i += 2;
        } else {
          noise_V[i] = noise_V[i+1] = rGIG_cpp(p_vec[i], a_vec[i], b_vec[i], rng());
          i++;
        }
      }
    }
  } else {
      // noise_V = rGIG_cpp(p_vec, a_vec, b_vec, rng());
      NoiseUtil::sample_V(noise_V, family, p_vec, a_vec, b_vec, rng);
  }
}

// posterior
Rcpp::List BlockModel::sampling(
  int n, int n_burnin, bool posterior, const SparseMatrix<double>& A
) {
  std::vector<VectorXd> AWs; // blockA * blockW
  std::vector<VectorXd> Ws; // return ZW (concatenated across latents)
  std::vector<VectorXd> cond_Ws; // return Z * E(W|Y,V)
  std::vector<VectorXd> Vs; // blockV
  std::vector<VectorXd> mn_Vs; // measurement nosie V
  rao_blackwell = true;

  if (!all_gaussian) burn_in(n_burnin);

  for (int i=0; i < n; i++) {
    if (posterior) {
      if (!all_gaussian) sample_cond_V(false);
      sampleW_VY();
      sample_cond_noise_V(true);
    } else {
      sample_uncond_V();
      sample_cond_noise_V(false);
    }

    // Collect A Z W for output
    if (n_latent > 0) {
      VectorXd AZW = VectorXd::Zero(n_obs);
      VectorXd ZW  = VectorXd::Zero(W_sizes);
      VectorXd ZcW = VectorXd::Zero(W_sizes);
      int woff = 0;
      for (int li = 0; li < n_latent; ++li) {
        const auto& Ai = latents[li]->getA();
        const auto& Zi = latents[li]->getZ();
        const auto& Wi = latents[li]->getW();
        const auto& cWi = latents[li]->get_cond_W();
        VectorXd ZiWi = Zi * Wi;
        VectorXd Zi_cWi = Zi * cWi;
        AZW.noalias() += Ai * ZiWi;
        ZW.segment(woff, ZiWi.size())  = ZiWi;
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

  return Rcpp::List::create(
    Rcpp::Named("AW") = AWs,
    Rcpp::Named("W") = Ws,
    Rcpp::Named("V") = Vs,
    Rcpp::Named("cond_W") = cond_Ws
  );
}

// fix parameter if converge
// void BlockModel::check_converge(vector<bool>& converge) {
//   int pos = 0;
//   for (std::vector<std::shared_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
//     int theta_len = (*it)->get_n_params();
//     vector<bool> sub_converge (converge.begin(), converge.begin() + theta_len);
//     (*it)->check_converge(sub_converge);
//     pos += theta_len;
//   }
// }

MatrixXd BlockModel::get_preconditioner(int strategy, double eps) {
  // Return cached preconditioner if available; do not compute here
  if (last_precond_valid && std::abs(eps - last_precond_eps) <= 1e-12 && strategy == last_precond_strategy) {
    return last_precond;
  }
  // Not computed/mismatch: identity-like fallback (caller should call compute_grad_and_hessian with with_precond=true)
  return VectorXd::Constant(n_params, 1e-5).asDiagonal();
}

// Main function for computing
void BlockModel::compute_grad_and_hessian(bool with_precond, double eps) {
if (debug) std::cout << "Start compute_grad_and_hessian"<< std::endl;
  auto t_total_start = std::chrono::steady_clock::now();
  long long t_sampleV_ms = 0, t_sampleW_ms = 0, t_rbtrace_ms = 0;
  long long t_build_s_ms = 0, t_set_s_ms = 0, t_dZ_ms = 0, t_grad_ms = 0;
  long long t_prec_latent_ms = 0, t_prec_ZGN_ms = 0, t_prec_merr_ms = 0;

  bool do_precond = with_precond;
  // Use the eps provided by the caller to keep caches consistent across layers
  double eps_use = eps;

  VectorXd latent_grad = VectorXd::Zero(n_la_params);
  VectorXd noise_grad  = VectorXd::Zero(n_params - n_la_params);
  MatrixXd precond_sum = MatrixXd::Zero(n_params, n_params);
  int      precond_count = 0;

  // RB: one pass using conditional W; Gibbs: n_gibbs passes using sampled W
  // Only collapse to a single RB pass for all-Gaussian models.
  // For non-Gaussian (e.g., NIG), we still need Gibbs over V even if RB is enabled.
  bool rb_all_gauss = (all_gaussian && rao_blackwell);
  int  n_pass = rb_all_gauss ? 1 : n_gibbs;
  // Use conditional mean of W for building observation score when RB is requested
  bool use_condW = rao_blackwell;

  if (rao_blackwell) assemble_dK();

  // For covariance of gradient samples (only meaningful when n_pass>1)
  MatrixXd grad_samples(n_params, std::max(1, n_pass));

  for (int i=0; i<n_pass; ++i) {
    if (rb_all_gauss) {
      if (i==0) {
        auto t_sw = std::chrono::steady_clock::now();
        sampleW_VY();     // compute QQ using current Z and V
        t_sampleW_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_sw).count();
        auto t_rb = std::chrono::steady_clock::now();
        compute_rb_trace(); // RB trace once (depends on QQ)
        t_rbtrace_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_rb).count();
      }
    } else {
      auto t_sv = std::chrono::steady_clock::now();
      // Avoid duplicating QQ factorization here; sampleW_VY() updates QQ
      sample_cond_V(false);
      t_sampleV_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_sv).count();
      auto t_sw = std::chrono::steady_clock::now();
      sampleW_VY(false);
      t_sampleW_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_sw).count();
      t_sv = std::chrono::steady_clock::now();
      sample_cond_noise_V();
      t_sampleV_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_sv).count();
      if (rao_blackwell) {
        auto t_rb = std::chrono::steady_clock::now();
        compute_rb_trace();
        t_rbtrace_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_rb).count();
      }
    }

    // Build observation score s = A^T D r for Z-chain
    auto t_bs = std::chrono::steady_clock::now();
    VectorXd s_full;
    if (!corr_measure && shared_sigma) {
      // D = diag(1 / sigma^2) when sigma is shared; use full residual
      VectorXd residual = get_residual(use_condW);
      VectorXd inv_sigma2 = noise_sigma.array().pow(-2).matrix();
      s_full = A.transpose() * inv_sigma2.asDiagonal() * residual;
    } else if (!corr_measure) {
      VectorXd residual = get_residual(use_condW);
      VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
      s_full = A.transpose() * inv_noise_SV.asDiagonal() * residual;
    } else {
      VectorXd residual = get_residual(use_condW);
      s_full = A.transpose() * Q_eps * residual;
    }
    t_build_s_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_bs).count();

    // Dispatch s to latents; compute dZ
    auto t_sets = std::chrono::steady_clock::now();
    int woff = 0;
    for (int li=0; li<n_latent; ++li) {
      int wi = latents[li]->get_W_size();
      latents[li]->set_obs_score(s_full.segment(woff, wi));
      woff += wi;
    }
    t_set_s_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_sets).count();

    // Aggregate gradients
    auto t_g = std::chrono::steady_clock::now();
    VectorXd current_grad = VectorXd::Zero(n_params);
    int pos = 0;
    for (auto it = latents.begin(); it != latents.end(); ++it) {
      int theta_len = (*it)->get_n_params();
      VectorXd gi = (*it)->get_grad(rao_blackwell);
      current_grad.segment(pos, theta_len) = gi;
      latent_grad .segment(pos, theta_len) += gi;
      pos += theta_len;
    }
    VectorXd noise_g = grad_theta_merr();
    current_grad.segment(n_la_params, n_merr) = noise_g;
    noise_grad.head(n_merr) += noise_g;
    if (!fix_flag[block_fix_beta]) {
      VectorXd beta_g = grad_beta();
      current_grad.segment(n_la_params + n_merr, n_feff) = beta_g;
      noise_grad.tail(n_feff) += beta_g;
    }
    t_grad_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_g).count();

    // Store gradient sample for covariance
    if (n_pass > 1) grad_samples.col(i) = current_grad;

    // Aggregate preconditioner blocks if requested
    if (do_precond) {
      // per-latent block preconditioners
      auto t_pl = std::chrono::steady_clock::now();
      int pos2 = 0;
      for (int li=0; li<n_latent; ++li) {
        int theta_len = latents[li]->get_n_params();
        MatrixXd Pi;
        try {
          Pi = latents[li]->preconditioner();
        } catch (...) {
          Pi = VectorXd::Constant(theta_len, latents[li]->get_V_size()).asDiagonal();
        }
        precond_sum.block(pos2, pos2, theta_len, theta_len) += Pi;
        pos2 += theta_len;
      }
      t_prec_latent_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_pl).count();

      // Z Gauss–Newton diagonal terms per latent (always diagonal wrt theta_K)
      auto t_pz = std::chrono::steady_clock::now();
      for (int li=0; li<n_latent; ++li) {
        SparseMatrix<double> ADA_i;
        const auto& Ai = latents[li]->getA();
        if (!corr_measure && shared_sigma) {
          ADA_i = Ai.transpose() * noise_V.cwiseInverse().asDiagonal() * Ai;
        } else if (!corr_measure) {
          VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
          ADA_i = Ai.transpose() * inv_noise_SV.asDiagonal() * Ai;
        } else {
          if (Q_eps.rows() == 0) {
            VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
            ADA_i = Ai.transpose() * inv_noise_SV.asDiagonal() * Ai;
          } else {
            ADA_i = Ai.transpose() * Q_eps * Ai;
          }
        }
        int base = 0; for (int k=0; k<li; ++k) base += latents[k]->get_n_params();
        VectorXd Wi_loc = use_condW ? latents[li]->get_cond_W() : latents[li]->getW();
        for (int j=0; j<latents[li]->get_n_theta_K(); ++j) {
          const auto& dZ_j = latents[li]->get_dZ(j);
          if (dZ_j.rows()==0 || dZ_j.cols()==0) continue;
          VectorXd g = dZ_j * Wi_loc;
          double hjj = g.transpose() * (ADA_i * g);
          precond_sum(base + j, base + j) += hjj;
        }
      }
      t_prec_ZGN_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_pz).count();

      // measurement/fixed-effects block
      auto t_pm = std::chrono::steady_clock::now();
      // 1) Analytic Hessian for measurement mu (noise mean) block:
      // H_mu = - B_mu^T diag(w) B_mu, with
      //   w_i = ((V_i - 1)^2) / V_i / sigma_i^2  (or / V_i if shared_sigma)
      if (n_theta_mu > 0 && !fix_flag[block_fix_theta_mu]) {
        VectorXd w;
        if (shared_sigma) {
          w = (noise_V.array() - 1.0).square().matrix().cwiseQuotient(noise_V);
        } else {
          VectorXd noise_SV = noise_sigma.array().square().matrix().cwiseProduct(noise_V);
          w = (noise_V.array() - 1.0).square().matrix().cwiseQuotient(noise_SV);
        }
        MatrixXd Hmu = (B_mu.transpose() * w.asDiagonal() * B_mu);
        // Place analytic H_mu at the top-left of the (merr+feff) block
        precond_sum.block(n_la_params + 0, n_la_params + 0, n_theta_mu, n_theta_mu) += Hmu;
      }

      // 1b) Analytic Hessian for measurement sigma and mu-sigma cross (when not shared_sigma)
      if (!shared_sigma && n_theta_sigma > 0 && !fix_flag[block_fix_theta_sigma]) {
        // e = residual = Y - mu' (V-1) - A Z W - X beta
        VectorXd e = get_residual(use_condW);
        // H_sigma = -2 B_sigma^T diag(e^2 / (sigma^2 V)) B_sigma
        VectorXd wsig = 2.0 * e.array().square().matrix()
                          .cwiseQuotient(noise_sigma.array().square().matrix().cwiseProduct(noise_V));
        MatrixXd Hsigma = (B_sigma.transpose() * wsig.asDiagonal() * B_sigma);
        // Place H_sigma just after the mu block within measurement corner
        precond_sum.block(n_la_params + n_theta_mu,
                          n_la_params + n_theta_mu,
                          n_theta_sigma,
                          n_theta_sigma) += Hsigma;
        // Cross H_{mu,sigma} = -2 B_mu^T diag(((V-1) ⊙ e) / (sigma^2 V)) B_sigma
        if (n_theta_mu > 0 && !fix_flag[block_fix_theta_mu]) {
          VectorXd wms = 2.0 * (noise_V - VectorXd::Ones(n_obs))
                               .cwiseProduct(e)
                               .cwiseQuotient(noise_sigma.array().square().matrix().cwiseProduct(noise_V));
          MatrixXd Hmu_sigma = (B_mu.transpose() * wms.asDiagonal() * B_sigma);
          precond_sum.block(n_la_params + 0,
                            n_la_params + n_theta_mu,
                            n_theta_mu,
                            n_theta_sigma) += Hmu_sigma;
          precond_sum.block(n_la_params + n_theta_mu,
                            n_la_params + 0,
                            n_theta_sigma,
                            n_theta_mu) += Hmu_sigma.transpose();
        }
      }

      // 1c) Analytic Hessian for measurement nu (no cross-terms)
      if (n_theta_nu > 0 && !fix_flag[block_fix_theta_nu]) {
        if (family != "normal") {
          // Use the same analytic form as latent: H_nu = - B_nu^T diag(nu ⊙ c) B_nu for NIG,
          // and the appropriate GAL/t variants handled inside NoiseUtil.
          MatrixXd Hnu = NoiseUtil::hess_theta_nu(family, B_nu, noise_nu, noise_V);
          precond_sum.block(n_la_params + n_theta_mu + n_theta_sigma,
                            n_la_params + n_theta_mu + n_theta_sigma,
                            n_theta_nu,
                            n_theta_nu) += Hnu;
        }
      }

      // Analytic cross-terms with beta
      // H_{mu,beta} = - B_mu^T diag((V-1)/(sigma^2 ∘ V)) X  (shared_sigma: drop sigma^2)
      // if (n_theta_mu > 0 && n_feff > 0 && !fix_flag[block_fix_theta_mu] && !fix_flag[block_fix_beta]) {
      //   VectorXd wmb;
      //   if (shared_sigma) {
      //     wmb = (noise_V.array() - 1.0).matrix().cwiseQuotient(noise_V);
      //   } else {
      //     VectorXd noise_SV = noise_sigma.array().square().matrix().cwiseProduct(noise_V);
      //     wmb = (noise_V.array() - 1.0).matrix().cwiseQuotient(noise_SV);
      //   }
      //   MatrixXd Hmu_beta = (B_mu.transpose() * wmb.asDiagonal() * X);
      //   // Place [mu,beta] and its transpose
      //   precond_sum.block(n_la_params + 0,
      //                     n_la_params + n_merr,
      //                     n_theta_mu,
      //                     n_feff) += Hmu_beta;
      //   precond_sum.block(n_la_params + n_merr,
      //                     n_la_params + 0,
      //                     n_feff,
      //                     n_theta_mu) += Hmu_beta.transpose();
      // }

      // H_{sigma,beta} = -2 B_sigma^T diag(e/(sigma^2 ∘ V)) X (only for non-shared sigma)
      // if (!shared_sigma && n_theta_sigma > 0 && n_feff > 0 && !fix_flag[block_fix_theta_sigma] && !fix_flag[block_fix_beta]) {
      //   VectorXd e_sb = get_residual(use_condW);
      //   VectorXd wsb = 2.0 * e_sb.cwiseQuotient(noise_sigma.array().square().matrix().cwiseProduct(noise_V));
      //   MatrixXd Hsigma_beta = (B_sigma.transpose() * wsb.asDiagonal() * X);
      //   precond_sum.block(n_la_params + n_theta_mu,
      //                     n_la_params + n_merr,
      //                     n_theta_sigma,
      //                     n_feff) += Hsigma_beta;
      //   precond_sum.block(n_la_params + n_merr,
      //                     n_la_params + n_theta_mu,
      //                     n_feff,
      //                     n_theta_sigma) += Hsigma_beta.transpose();
      // }

      // 1d) Analytic Hessian for beta (fixed effects): H_beta = - X^T Sigma^{-1} X
      if (n_feff > 0 && !fix_flag[block_fix_beta]) {
        MatrixXd Hbeta;
        // if (shared_sigma) {
        //   // Sigma^{-1} = I in the transformed residual form used by grad_beta()
        //   Hbeta = (X.transpose() * X);
        // } else {
        //   VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
        //   Hbeta = (X.transpose() * inv_noise_SV.asDiagonal() * X);
        // }
        Hbeta = MatrixXd::Identity(n_feff, n_feff) * n_obs;
        precond_sum.block(n_la_params + n_merr,
                          n_la_params + n_merr,
                          n_feff,
                          n_feff) += Hbeta;
      }

      t_prec_merr_ms += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_pm).count();

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
  noise_grad  /= std::max(1, n_pass);

  VectorXd avg_gradient = VectorXd::Zero(n_params);
  avg_gradient.head(n_la_params) = latent_grad;
  avg_gradient.tail(n_params - n_la_params) = noise_grad;

  // Finalize cached preconditioner
  if (do_precond) {
    if (precond_count > 0) {
      last_precond = (1.0 / precond_count) * precond_sum;
      last_precond += VectorXd::Constant(n_params, 1e-5).asDiagonal();
      last_precond_valid = true;
      last_precond_eps = eps_use;
      last_precond_strategy = precond_strategy_;
    } else {
      last_precond_valid = false;
    }
  }

  last_gradient = avg_gradient;
  last_grad_valid = true;

  if (debug) {
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_total_start).count();
    std::cout << "[block] compute_grad_and_hessian timing (ms): total=" << total_ms
              << ", sampleV=" << t_sampleV_ms
              << ", sampleW=" << t_sampleW_ms
              << ", rb_trace=" << t_rbtrace_ms
              << ", build_s=" << t_build_s_ms
              << ", set_s=" << t_set_s_ms
              << ", dZ=" << t_dZ_ms
              << ", grad=" << t_grad_ms
              << ", prec_latent=" << t_prec_latent_ms
              << ", prec_ZGN=" << t_prec_ZGN_ms
              << ", prec_merr=" << t_prec_merr_ms
              << std::endl;
  }
}
  

void BlockModel::compute_rb_trace() {
if (debug) std::cout << "start compute trace" << std::endl;
  auto t_start = std::chrono::steady_clock::now();
  int n = 0;
  int woff = 0; // offset in W-space for embedding latent-local pieces
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
  for (int i=0; i < n_latent; i++) {
    VectorXd rb_trace_K (latents[i]->get_n_theta_K());
    VectorXd rb_trace_sigma (latents[i]->get_n_theta_sigma());

    // compute for K: tr(QQ^-1 dK^T diag(1/SV) K)
    for (int j=0; j < latents[i]->get_n_theta_K(); j++) {
      SparseMatrix<double> T = block_dK[i][j].transpose() * inv_SV.asDiagonal() * K;
      rb_trace_K[j] = chol_QQ.trace(T, rng());
    }

    // compute for sigma: tr(Q^-1 K B_sigma.col(j)/SV K^T) for non-fixed theta_sigma
    vector<bool> fix_theta_sigma_vec = latents[i]->get_theta_unfixed_sigma();
    int pos=0;
    for (int j=0; j < latents[i]->get_n_theta_sigma(); j++) {
      if (fix_theta_sigma_vec[pos]) pos += 1; // find non-fixed theta_sigma position

      // build B_sigma_col_j (consider all latents)
      VectorXd BSigma_col_over_SV = VectorXd::Zero(V_sizes);
      BSigma_col_over_SV.segment(n, latents[i]->get_V_size()) = latents[i]->get_BSigma_col(pos);
      BSigma_col_over_SV = BSigma_col_over_SV.cwiseProduct(inv_SV);

      SparseMatrix<double> T = K.transpose() * BSigma_col_over_SV.asDiagonal() * K;
      rb_trace_sigma[j] = chol_QQ.trace(T, rng());
      pos += 1;
    }

    // Add Z-related RB trace: T = (dZ_j)^T A_i^T D A_i Z_i
    // where D is measurement precision (depends on noise settings)
    // Matches: T = dZ^T A^T D A Z
    // where D is measurement precision (depends on noise settings)
    // Precompute A_i^T D A_i once per latent
    SparseMatrix<double> ADA_i;
    {
      const auto& Ai = latents[i]->getA();
      if (!corr_measure && shared_sigma) {
        // D = diag(1 / V)
        ADA_i = Ai.transpose() * noise_V.cwiseInverse().asDiagonal() * Ai;
      } else if (!corr_measure) {
        // D = diag(1 / (sigma^2 V))
        VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
        ADA_i = Ai.transpose() * inv_noise_SV.asDiagonal() * Ai;
      } else {
        // Correlated case: D = Q_eps
        ADA_i = Ai.transpose() * Q_eps * Ai;
      }
    }

    // Compute and accumulate per-parameter Z traces
    for (int j=0; j<latents[i]->get_n_theta_K(); ++j) {
      const auto& dZ_j = latents[i]->get_dZ(j);
      if (dZ_j.rows()==0 || dZ_j.cols()==0) continue; // skip if operator doesn't provide Z
      const auto& Zi = latents[i]->getZ();
      // Local block in latent i's W-space
      SparseMatrix<double> Tloc = dZ_j.transpose() * ADA_i * Zi;
      // Embed into full W-space
      SparseMatrix<double> T(W_sizes, W_sizes);
      setSparseBlock(&T, woff, woff, Tloc);
      // Accumulate into rb_trace_K (same sign as K-term; latent subtracts rb_trace_K later)
      rb_trace_K[j] += chol_QQ.trace(T, rng());
    }

    latents[i]->set_rb_trace(rb_trace_K, rb_trace_sigma);
    n += latents[i]->get_V_size();
    woff += latents[i]->get_W_size();
  }

  // compute for theta_sigma
  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());
  for (int j=0; j < n_theta_sigma; j++) {
    // Include Z: T = (A Z)^T diag(B_sigma_j / noise_SV) (A Z)
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    SparseMatrix<double> T = AZ.transpose() * B_sigma.col(j).cwiseQuotient(noise_SV).asDiagonal() * AZ;
    rb_trace_noise_sigma[j] = chol_QQ.trace(T, rng());
  }
  if (debug) {
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_start).count();
    std::cout << "after compute trace (" << ms << " ms)" << std::endl;
  }
}

void BlockModel::assemble_dK() {
  int nrow = 0; int ncol = 0;
  for (int i=0; i < n_latent; i++) {
      for (int j=0; j < latents[i]->get_n_theta_K(); j++) {
          setSparseBlock(&block_dK[i][j], nrow, ncol, latents[i]->get_dK(j));
      }
      nrow += latents[i]->get_V_size();
      ncol += latents[i]->get_W_size();
  }
}

// generate output to R
Rcpp::List BlockModel::output() const {
  Rcpp::List latents_output;
  for (int i=0; i < n_latent; i++) {
    latents_output.push_back((*latents[i]).output());
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("noise") = Rcpp::List::create(
      Rcpp::Named("noise_type")   = family,
      Rcpp::Named("theta_mu")     = theta_mu,
      Rcpp::Named("theta_sigma")  = theta_sigma,
      Rcpp::Named("theta_nu")     = theta_nu,
      Rcpp::Named("V")            = noise_V,
      Rcpp::Named("rho")          = rho
    ),
    Rcpp::Named("feff")            = beta,
    // Rcpp::Named("sampling_time")   = sampling_time.count(),
    Rcpp::Named("models")          = latents_output
    // Rcpp::Named("log_likelihood")  = all_gaussian ? -log_likelihood() : 0
  );

  return out;
}


SparseMatrix<double> BlockModel::get_sqrt_AtSVA() const {
  VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);  
  if (!corr_measure && shared_sigma) {
    // H = sqrt(Sigma^{-1}) * A Z
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    return noise_V.cwiseSqrt().cwiseInverse().asDiagonal() * AZ;
  } else if (!corr_measure) {
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    return inv_noise_SV.cwiseSqrt().asDiagonal() * AZ;
  } else {
    SparseMatrix<double> AZ(n_obs, W_sizes);
    int col = 0;
    for (int li=0; li<n_latent; ++li) {
      SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
      setSparseBlock(&AZ, 0, col, AiZi);
      col += latents[li]->get_W_size();
    }
    return sqrt_Rinv * inv_noise_SV.cwiseSqrt().asDiagonal() * AZ;
  }
}

// update Q_eps, dQ_eps, and compute trace as sum_ij dQ_ij * Q^-1_ij
void BlockModel::update_Q_eps(double rho) {
    // update Q_eps
    for (int i=0; i < Q_eps.outerSize(); i++) {
      for (SparseMatrix<double>::InnerIterator it(Q_eps, i); it; ++it) {
        if (it.row() == it.col()) {
          int idx = it.row();
          it.valueRef() = 1.0/(pow(noise_sigma(idx), 2) * noise_V(idx));
          if (has_correlation[idx]) it.valueRef() /= (1-rho*rho);
        } else {
          double tmp = noise_sigma(it.row()) * noise_sigma(it.col()) * sqrt(noise_V(it.row()) * noise_V(it.col()));
          it.valueRef() = -rho / ((1-rho*rho) * tmp);
        }
      }
    }

    // update sqrt_Rinv
    for (int i=0; i < sqrt_Rinv.outerSize(); i++) {
      for (SparseMatrix<double>::InnerIterator it(sqrt_Rinv, i); it; ++it) {
        if (it.row() == it.col()) {
          int idx = it.row();
          it.valueRef() = has_correlation[idx] ? 
            (sqrt(1+rho)+sqrt(1-rho))/(2*sqrt(1-rho*rho)) : 1;
        } else {
          it.valueRef() = (-sqrt(1+rho)+sqrt(1-rho))/(2*sqrt(1-rho*rho));
        }
      }
    }

    // update dQ_eps, and compute trace as sum_ij dQ_ij * Q^-1_ij
    for (int i=0; i < dQ_eps.outerSize(); i++) {
      for (SparseMatrix<double>::InnerIterator it(dQ_eps, i); it; ++it) {
        if (it.row() == it.col()) {
          int idx = it.row();
          double tmp = pow((1-rho*rho) * noise_sigma(idx), 2) * noise_V(idx);
          it.valueRef() = 2.0*rho / tmp;
        } else {
          int r = it.row(); int c = it.col();
          double tmp = pow(1-rho*rho, 2) * noise_sigma(r) * noise_sigma(c) * sqrt(noise_V(r) * noise_V(c));
          it.valueRef() = -(1+rho*rho) / tmp;
        }
      }
    }

    // compute logdet of Q_eps
    // lhs = Q_eps.logdet();
    VectorXd SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_prevV);
    double lhs=0; int i=0;
    while (i < n_obs) {
      if (has_correlation[i]) {
        lhs += - log((1-rho*rho) * SV[i] * SV[i+1]);
        i+=2;
      } else {
        lhs += - log(SV[i]);
        i+=1;
      }
    }
    VectorXd res = get_residual();
    double rhs = res.dot(Q_eps * res);
}

void BlockModel::update_QQ() {
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());

  // update Q and QQ
  Q = K.transpose() * inv_SV.asDiagonal() * K;
  // Build H = sqrt(D) * A Z consistent with get_sqrt_AtSVA(), then QQ = Q + H^T H
  SparseMatrix<double> AZ(n_obs, W_sizes);
  int col = 0;
  for (int li=0; li<n_latent; ++li) {
    SparseMatrix<double> AiZi = latents[li]->getA() * latents[li]->getZ();
    setSparseBlock(&AZ, 0, col, AiZi);
    col += latents[li]->get_W_size();
  }
  if (!corr_measure && shared_sigma) {
    // D = diag(1 / V)
    SparseMatrix<double> H = noise_V.cwiseSqrt().cwiseInverse().asDiagonal() * AZ;
    QQ = Q + H.transpose() * H;
  } else if (!corr_measure) {
    // D = diag(1 / (sigma^2 V))
    VectorXd inv_noise_SV = noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V);
    SparseMatrix<double> H = inv_noise_SV.cwiseSqrt().asDiagonal() * AZ;
    QQ = Q + H.transpose() * H;
  } else {
    // Correlated case: D = Q_eps
    QQ = Q + AZ.transpose() * Q_eps * AZ;
  }
  chol_QQ.compute(QQ);
}
