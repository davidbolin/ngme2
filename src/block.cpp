// Implementation for block model and block_rep

#include "block.h"
#include "prior.h"
#include "sample_rGIG.h"
#include <random>
#include <cmath>
#include <iterator>

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
  // dK            (V_sizes, W_sizes)
  // d2K           (V_sizes, W_sizes)
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
    latent_in["n_trace_iter"] = n_trace_iter;
    unsigned long latent_seed = rng();
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
    fix_flag[block_fix_theta_sigma]  = Rcpp::as<bool> (noise_in["fix_theta_sigma"]);
    fix_flag[blcok_fix_V]            = Rcpp::as<bool> (noise_in["fix_V"]);
    fix_flag[block_fix_nu]           = Rcpp::as<bool> (noise_in["fix_nu"]);

    B_mu          = (Rcpp::as<MatrixXd>      (noise_in["B_mu"])),
    theta_mu      = (Rcpp::as<VectorXd>      (noise_in["theta_mu"])),
    n_theta_mu    = (theta_mu.size()),

    B_sigma       = (Rcpp::as<MatrixXd>      (noise_in["B_sigma"])),
    theta_sigma   = (Rcpp::as<VectorXd>      (noise_in["theta_sigma"])),
    n_theta_sigma = (theta_sigma.size()),

    B_nu          = (Rcpp::as<MatrixXd>      (noise_in["B_nu"])),
    theta_nu      = (Rcpp::as<VectorXd>      (noise_in["theta_nu"])),
    n_theta_nu    = (theta_nu.size()),

    rb_trace_noise_sigma = VectorXd::Zero(n_theta_sigma),

    family = Rcpp::as<string>  (noise_in["noise_type"]);
    noise_mu = B_mu * theta_mu;
    noise_sigma = (B_sigma * theta_sigma).array().exp();
    noise_nu = (B_nu * theta_nu).array().exp();

    rho = Rcpp::as<VectorXd> (noise_in["rho"]); n_rho = rho.size();
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

      SparseMatrix<double> dQ_lower (n_obs, n_obs);
      dQ_lower.setFromTriplets(dQ_eps_triplet.begin(), dQ_eps_triplet.end());
      dQ_eps = dQ_lower.selfadjointView<Lower>();

    // Q_eps_solver.analyzePattern(Q_eps);
// std::cout << "Q_eps: \n" << Q_eps << std::endl;
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
    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    if (!corr_measure) {
      QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * A;
    }
    else{
      QQ = Q + A.transpose() * Q_eps * A;
    }

    // set N for trace estimator
    chol_QQ.set_N(n_trace_iter);

    chol_Q.analyze(Q);
    chol_QQ.analyze(QQ);
    LU_K.analyzePattern(K);
  }
if (debug) std::cout << "After init solver " << std::endl;

  // 8. optimizer related
  // stepsizes = VectorXd::Constant(n_params, stepsize);
  steps_to_threshold = VectorXd::Constant(n_params, 0);
  indicate_threshold = VectorXd::Constant(n_params, 0);

if (debug) std::cout << "After init solver && before sampleW_V" << std::endl;

  if (n_latent > 0 && init_sample_W) {
    sampleW_V();
    sampleW_V();
  }

  // if (all_gaussian) rao_blackwell = true;
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
if (debug) std::cout << "End Block Constructor" << std::endl;
}

  void BlockModel::burn_in(int iterations) {
      // sample_uncond_V();
      // sample_uncond_noise_V();
      for (int i=0; i < iterations; i++) {
// std::cout << "burn in iteration i= " << i  << std::endl;
          sample_cond_V();
// std::cout << "sample cond V done" << std::endl;
          sampleW_VY(true);
// std::cout << "sample W done" << std::endl;
          sample_cond_noise_V();
// std::cout << "sample cond_noise_V done" << std::endl;
      }
// std::cout << "burn in done "  << std::endl;
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
  // cond_W = W;

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
// if (debug) std::cout << "starting sampling W." << std::endl;
// double time = 0;
  if (n_latent==0) return;

  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());

  // init Q and QQ, Q is updated every iteration, so does QQ
  Q = K.transpose() * inv_SV.asDiagonal() * K;
  VectorXd residual_part = get_residual_part();
  VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean();

  if (!corr_measure) {
    QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * A;
    M += A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * residual_part;
  }
  else{
    QQ = Q + A.transpose() * Q_eps * A;
    M += A.transpose() * Q_eps * residual_part;
  }
// std::cout << "Q: \n" << Q << std::endl;
// std::cout << "QQ: \n" << QQ << std::endl;

// std::cout << "M = " << M << std::endl;
  VectorXd z = NoiseUtil::rnorm_vec(W_sizes, 0, 1, rng());

  // sample W ~ N(QQ^-1*M, QQ^-1)
  chol_QQ.compute(QQ);
  VectorXd W = chol_QQ.rMVN(M, z);
  setW(W);

  if (rao_blackwell && !burn_in) {
    // compute E(W|V,Y) i.e. QQ^-1 M
    W = chol_QQ.solve(M);
// std::cout << "cond W = " << W.transpose() << std::endl;
    set_cond_W(W);
  }

// std::cout << "size of W and time of sampling is " << W.size() << " " << time << std::endl;

if (debug) std::cout << "Finish sampling W" << std::endl;
}

// ---------------- get, set update gradient ------------------
// order is Latent, merr, feff
VectorXd BlockModel::get_parameter() {
if (debug) std::cout << "Start get_parameter"<< std::endl;
    VectorXd thetas (n_params);
    int pos = 0;
    for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      VectorXd theta = (*it)->get_parameter();
      thetas.segment(pos, theta.size()) = theta;
      pos += theta.size();
    }
    thetas.segment(n_la_params, n_merr) = get_theta_merr();

    if (!fix_flag[block_fix_beta] ) {
      thetas.segment(n_la_params + n_merr, n_feff) = beta;
    }

if (debug) std::cout << "Finish get_parameter"<< std::endl;
    return thetas;
}

// avg over gibbs samples
VectorXd BlockModel::grad() {
if (debug) std::cout << "Start block gradient"<< std::endl;
long long time_compute_g = 0;
long long time_sample_w = 0;
// auto timer_computeg = std::chrono::steady_clock::now();
// auto timer_sampleW = std::chrono::steady_clock::now();
// time_sample_w += since(timer_sampleW).count();

  VectorXd latent_grad = VectorXd::Zero(n_la_params);
  VectorXd noise_grad = VectorXd::Zero(n_params-n_la_params);

  if (rao_blackwell) assemble_dK(); // for computing trace
  if (all_gaussian && rao_blackwell) {
    // compute RB version of gradient
    sampleW_VY(); // QQ.compute
    compute_rb_trace();
    // get the grad (Latent, merr, feff)
    int pos = 0;
    for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      int theta_len = (*it)->get_n_params();
      latent_grad.segment(pos, theta_len) += (*it)->get_grad(rao_blackwell);
      pos += theta_len;
    }
    noise_grad.head(n_merr) += grad_theta_merr();
    if (!fix_flag[block_fix_beta]) {
      noise_grad.tail(n_feff) += grad_beta();
    }
  } else { // Running Gibbs sampling
    for (int i=0; i < n_gibbs; i++) {
// std::chrono::steady_clock::time_point startTime, endTime; startTime = std::chrono::steady_clock::now();
      // bool QQ_update = (i==0) || family != "normal";
      sample_cond_V();
      sampleW_VY(false);
      sample_cond_noise_V();
      if (rao_blackwell) compute_rb_trace();

      int pos = 0;
      for (std::vector<std::shared_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
        int theta_len = (*it)->get_n_params();
        latent_grad.segment(pos, theta_len) += (*it)->get_grad(rao_blackwell);
        pos += theta_len;
      }

// endTime = std::chrono::steady_clock::now(); sampling_time += std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);std::cout << "!!sampling time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << std::endl;
      noise_grad.head(n_merr) += grad_theta_merr();
      if (!fix_flag[block_fix_beta]) {
        noise_grad.tail(n_feff) += grad_beta();
      }
    }
    noise_grad  = (1.0/n_gibbs) * noise_grad;
    latent_grad = (1.0/n_gibbs) * latent_grad;
  }

  VectorXd avg_gradient = VectorXd::Zero(n_params);
  avg_gradient.head(n_la_params) = latent_grad;
  avg_gradient.tail(n_params-n_la_params) = noise_grad;

  // EXAMINE the gradient to change the stepsize
  if (reduce_var) examine_gradient();

if (debug) {
  std::cout << "avg time for compute grad (ms): " << time_compute_g / n_gibbs << std::endl;
  std::cout << "avg time for sampling W(ms): " << time_sample_w / n_gibbs << std::endl;
  std::cout << "Finish block gradient"<< std::endl;
}

  return avg_gradient;
}

void BlockModel::set_parameter(const VectorXd& Theta) {
if (debug) std::cout << "Start set_parameter"<< std::endl;
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
if (debug) std::cout << "Finish set_parameter"<< std::endl;
}

// sample W|V
void BlockModel::sampleW_V()
{
// if (debug) std::cout << "starting sampling W." << std::endl;
  if(n_latent==0) return;
  std::normal_distribution<double> rnorm {0,1};

  // sample KW ~ N(mu*(V-h), diag(V))
  VectorXd SV = getSV();
  Eigen::VectorXd KW = VectorXd::Zero(V_sizes);
  for (int i=0; i < V_sizes; i++) {
    // KW[i] = R::rnorm(0, sqrt(SV[i]));
    KW[i] = rnorm(rng) * sqrt(SV[i]);
  }
  KW = getMean() + KW;

  VectorXd W = VectorXd::Zero(W_sizes);
  if (V_sizes == W_sizes) {
    LU_K.factorize(K);
    int success = LU_K.info();
    if (success == 0) W = LU_K.solve(KW);
  } else {
    // SparseMatrix<double> Q = K.transpose() * K;
    // chol_Q.compute(Q);
    // W = chol_Q.solve(K.transpose() * KW);
    // W = K.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(KW);
  }
// std::cout << "K = " << K << std::endl;
  setW(W);
}

// --------- Fiexed effects and Measurement Error ---------------
VectorXd BlockModel::grad_beta() {
  VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());

  VectorXd residual = get_residual(rao_blackwell); // + X * beta;
  VectorXd grads = X.transpose() * noise_inv_SV.asDiagonal() * residual;
  //  * residual.cwiseQuotient(noise_sigma);

  // MatrixXd hess = X.transpose() * noise_inv_SV.asDiagonal() * X;
  // grads = grads / A.rows();
  // grads = hess.ldlt().solve(grads);

  // std::cout << "grads of beta=" << grads << std::endl;
    return -grads;
}

VectorXd BlockModel::grad_theta_mu() {
  // VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());
  // MatrixXd noise_X = (-VectorXd::Ones(n_obs) + noise_V).asDiagonal() * B_mu;

  // VectorXd grad = noise_X.transpose() * noise_inv_SV.asDiagonal() * residual;
  // MatrixXd hess = noise_X.transpose() * noise_inv_SV.asDiagonal() * noise_X;
  // grad = hess.ldlt().solve(grad);

  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());

  // VectorXd residual = get_residual(false);
  VectorXd residual = get_residual(rao_blackwell);
  VectorXd grad = VectorXd::Zero(n_theta_mu);
  for (int l=0; l < n_theta_mu; l++) {

      // LU_K.factorize(getK());
      // VectorXd tmp = residual + LU.solve(-h + )
      grad(l) = (noise_V - VectorXd::Ones(n_obs)).cwiseProduct(B_mu.col(l).cwiseQuotient(noise_SV)).dot(residual);

      // add prior
      // grad(l) += PriorUtil::d_log_dens(prior_mu_type, prior_mu_param, theta_mu(l));
  }
  return -grad;
}

VectorXd BlockModel::grad_theta_sigma() {
  VectorXd grad = VectorXd::Zero(n_theta_sigma);
  VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);
  // grad = B_sigma.transpose() * (-0.5 * VectorXd::Ones(n_obs) + residual.array().pow(2).matrix().cwiseQuotient(noise_SV));
  // std::cout << "get_cond_W = " << get_cond_W().mean() << std::endl;
  // std::cout << "get_W = " << getW().mean() << std::endl;
  // std::cout << "trace = " << rb_trace_noise_sigma << std::endl;
  // std::cout << "-----  " << std::endl;

  // VectorXd residual = get_residual(false);
  VectorXd residual = get_residual(rao_blackwell);
  VectorXd vsq = (residual).array().pow(2).matrix().cwiseQuotient(noise_SV);
  VectorXd tmp1 = vsq - VectorXd::Ones(n_obs);
  grad = B_sigma.transpose() * tmp1;
  if (rao_blackwell) grad += rb_trace_noise_sigma;

  // add prior
  // for (int l=0; l < n_theta_sigma; l++) {
  //     grad(l) += PriorUtil::d_log_dens(prior_sigma_type, prior_sigma_param, theta_sigma(l));
  // }

  return -grad;
}

VectorXd BlockModel::get_theta_merr() const {
  VectorXd theta_merr = VectorXd::Zero(n_merr);

  theta_merr.segment(0, n_theta_mu) = theta_mu;
  theta_merr.segment(n_theta_mu, n_theta_sigma) = theta_sigma;
  if (family != "normal")
    theta_merr.tail(n_theta_nu) = theta_nu;

// std::cout << " rho === " << rho << std::endl;
  if (corr_measure) theta_merr(n_merr-1) = rho2th(rho(0));

  return theta_merr;
}

VectorXd BlockModel::grad_theta_merr() {
  VectorXd grad = VectorXd::Zero(n_merr);

  if (!fix_flag[block_fix_theta_mu])     grad.segment(0, n_theta_mu) = grad_theta_mu();
  if (!fix_flag[block_fix_theta_sigma])  grad.segment(n_theta_mu, n_theta_sigma) = grad_theta_sigma();
  if (!fix_flag[block_fix_nu] && family != "normal") {
    grad.segment(n_theta_mu + n_theta_sigma, n_theta_nu) = NoiseUtil::grad_theta_nu(family, B_nu, noise_nu, noise_V, noise_prevV);

    // add prior
    // grad(n_theta_mu + n_theta_sigma) -= PriorUtil::d_log_dens(prior_nu_type, prior_nu_param, noise_nu);
  }

  // grad of theta_rho
  if (corr_measure) {
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
  theta_mu = theta_merr.segment(0, n_theta_mu);
  theta_sigma = theta_merr.segment(n_theta_mu, n_theta_sigma);
  if (family != "normal") {
    theta_nu = theta_merr.segment(n_theta_mu + n_theta_sigma, n_theta_nu);
    if (theta_nu(0) > log(1e4)) theta_nu(0) = 1e4;
  }

  // update mu, sigma
  noise_mu    = (B_mu * theta_mu);
  noise_sigma = (B_sigma * theta_sigma).array().exp();
  noise_nu    = (B_nu * theta_nu).array().exp();

  // update rho, and Q_eps
  if (corr_measure) {
    rho(0) = th2rho(theta_merr(n_merr-1));
    // update Q_eps
    for (int i=0; i < Q_eps.outerSize(); i++) {
      for (SparseMatrix<double>::InnerIterator it(Q_eps, i); it; ++it) {
        if (it.row() == it.col()) {
          int idx = it.row();
          it.valueRef() = 1.0/(pow(noise_sigma(idx), 2) * noise_V(idx));
          if (has_correlation[idx]) it.valueRef() /= (1-rho(0)*rho(0));
        } else {
          double tmp = noise_sigma(it.row()) * noise_sigma(it.col()) * sqrt(noise_V(it.row()) * noise_V(it.col()));
          it.valueRef() = -rho(0) / ((1-rho(0)*rho(0)) * tmp);
        }
      }
    }

    // update dQ_eps, and compute trace as sum_ij dQ_ij * Q^-1_ij
    for (int i=0; i < dQ_eps.outerSize(); i++) {
      for (SparseMatrix<double>::InnerIterator it(dQ_eps, i); it; ++it) {
        if (it.row() == it.col()) {
          int idx = it.row();
          double tmp = pow((1-rho(0)*rho(0)) * noise_sigma(idx), 2) * noise_V(idx);
          it.valueRef() = 2.0*rho(0) / tmp;
        } else {
          int r = it.row(); int c = it.col();
          double tmp = pow(1-rho(0)*rho(0), 2) * noise_sigma(r) * noise_sigma(c) * sqrt(noise_V(r) * noise_V(c));
          it.valueRef() = -(1+rho(0)*rho(0)) / tmp;
        }
      }
    }

  }
// show the construction
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
  std::vector<VectorXd> Ws; // blockW
  std::vector<VectorXd> Vs; // blockV
  std::vector<VectorXd> mn_Vs; // measurement nosie V

  burn_in(n_burnin);

  for (int i=0; i < n; i++) {
    if (posterior) {
      sample_cond_V();
      sampleW_VY();
      sample_cond_noise_V(true);
    } else {
      sample_uncond_V();
      sampleW_V();
      sample_cond_noise_V(false);
    }

    AWs.push_back(A * getW());
    Ws.push_back(getW());
    Vs.push_back(getV());
    // mn_Vs.push_back(var.getV());
  }

  return Rcpp::List::create(
    Rcpp::Named("AW") = AWs,
    Rcpp::Named("W") = Ws,
    Rcpp::Named("V") = Vs
    // Rcpp::Named("noise_V") = mn_Vs
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


// provide stepsize
void BlockModel::examine_gradient() {

    // examine if the gradient under the threshold
    // for (int i=0; i < n_params; i++) {
    //     if (abs(gradients(i)) < threshold) {
    //         indicate_threshold(i) = 1;
    //     }      // mark if under threshold
    //     if (!indicate_threshold(i)) steps_to_threshold(i) = counting; // counting up
    // }

    // counting += 1;
    // stepsizes = (VectorXd::Constant(n_params, counting) - steps_to_threshold).cwiseInverse().array().pow(reduce_power);

    // // finish opt fo latents
    // for (int i=0; i < n_latent; i++) {
    //     for (int j=0; j < latent_para; j++) {
    //         int index = latent_para*i + j;

    //         // if (counting - steps_to_threshold(index) > 100)
    //         if (abs(gradients(index)) < termination)
    //             latents[i]->finishOpt(j);
    //     }
    // }

// if (debug) {
//     std::cout << "steps=" << steps_to_threshold <<std::endl;
//     std::cout << "gradients=" << gradients <<std::endl;
//     std::cout << "stepsizes=" << stepsizes <<std::endl;
// }
}

MatrixXd BlockModel::precond(int strategy, double eps) {
if (debug) std::cout << "start precond"<< std::endl;
  MatrixXd precond = MatrixXd::Zero(n_params, n_params);
  // 1. Preconditioner for Latents
  int index_params = 0;
  for (int i=0; i < n_latent; i++) {
    int n_la = latents[i]->get_n_params();
    if (strategy == 0) {
      // No preconditioner
      precond.block(index_params, index_params, n_la, n_la) =
        VectorXd::Constant(n_la, latents[i]->get_V_size()).asDiagonal();
        // VectorXd::Constant(n_la, n_obs).asDiagonal();
    } else if (strategy == 1) {
      // Fast preconditioner
      precond.block(index_params, index_params, n_la, n_la) = latents[i]->precond(false, eps);
    } else if (strategy == 2) {
      // Full preconditioner
      precond.block(index_params, index_params, n_la, n_la) = latents[i]->precond(true, eps);
    }
    index_params += n_la;
  }
if (debug) std::cout << "after latents precond"<< std::endl;

  // 2. Preconditioner for fixed effects and measurement error
  if (strategy == 0) {
    // No preconditioner (or the 1st iteration)
    precond.bottomRightCorner(n_merr + n_feff, n_merr + n_feff) = VectorXd::Constant(n_merr + n_feff, n_obs).asDiagonal();
  } else {
    // have preconditioner
    VectorXd v (n_merr + n_feff);
    VectorXd th_rho = rho.unaryExpr(std::ref(rho2th));
    v << theta_mu, theta_sigma, theta_nu, th_rho, beta;

    MatrixXd hess_merr_feff = num_h_no_latent(v, eps);
    precond.bottomRightCorner(n_merr + n_feff, n_merr + n_feff) = hess_merr_feff;
// std::cout << "merr feff precond = " << hess_merr_feff <<std::endl;
  }
  // add small eps to diagonal
  precond += VectorXd::Constant(n_params, 1e-5).asDiagonal();

// std::cout << "block precond =" << precond <<std::endl;
  return precond;
}

// precond_fast: numerical hessian for c(theta_mu, theta_sigma, nu, rho, feff)
MatrixXd BlockModel::num_h_no_latent(const VectorXd& v, double eps) {
	int n = v.size();
	MatrixXd hessian(n, n);
	double original_val = logd_no_latent(v);
	VectorXd f_v (n);
	// compute f_v = logd_no_latent( + eps * e_i)
	for (int i=0; i < n; i++) {
		VectorXd tmp_v = v; tmp_v(i) += eps;
		f_v(i) = logd_no_latent(tmp_v);
	}

	// compute H_ij = d2 f / dxi dxj
	for (int i=0; i < n; i++) {
		for (int j=0; j <= i; j++) {
			VectorXd tmp_vij = v; tmp_vij(i) += eps; tmp_vij(j) += eps;
			double f_vij = logd_no_latent(tmp_vij);
			double h_ij = (f_vij - f_v(i) - f_v(j) + original_val) / (eps * eps);
      hessian(i, j) = h_ij;
		}
	}

	// fill in the lower triangular part
	for (int i=0; i < n; i++) {
		for (int j=0; j < i; j++) {
			hessian(j, i) = hessian(i, j);
		}
	}
	return hessian;
}

// v = mu, sigma, nu, rho, feff
// residual ~ N(-mu + mu V, .. )
double BlockModel::logd_no_latent(const VectorXd& v) {
	VectorXd theta_mu    = v.segment(0, n_theta_mu);
	VectorXd theta_sigma = v.segment(n_theta_mu, n_theta_sigma);
	VectorXd theta_nu    = v.segment(n_theta_mu + n_theta_sigma, n_theta_nu);
  VectorXd th_rho = v.segment(n_theta_mu + n_theta_sigma + n_theta_nu, n_rho);
  VectorXd beta = v.tail(n_feff);

  // compute mu sigma, nu
  VectorXd noise_mu    = (B_mu * theta_mu);
  VectorXd noise_sigma = (B_sigma * theta_sigma).array().exp();

  // pi(Y|W, V) (no correction for noise!)
  // residual but using prevV
	VectorXd tmp = Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + noise_prevV).cwiseProduct(noise_mu);
  VectorXd SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_prevV);

  double logd_res = 0;
  if (corr_measure) {
    // No need to copy construct
    // SparseMatrix<double> Q_eps_tmp = Q_eps;
    // update rho and Q_eps
    double rho = th_rho.unaryExpr(std::ref(th2rho))(0);
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
    // compute logdet of Q_eps
    // lhs = Q_eps.logdet();
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
    double rhs = tmp.dot(Q_eps * tmp);
    logd_res = 0.5 * (lhs - rhs);
  } else {
    double lhs = SV.array().log().sum();
    double rhs = tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    logd_res = - 0.5 * (lhs + rhs);
  }

  // pi(V)
  double logd_V = 0;
  VectorXd h = VectorXd::Ones(n_obs);
  if (family != "normal") {
    logd_V = NoiseUtil::log_density(family, noise_V, h, B_nu, theta_nu, FALSE);
  }

// std::cout << "logd_res=" << logd_res << std::endl;
  return  -(logd_res + logd_V);
}

void BlockModel::compute_rb_trace() {
if (debug) std::cout << "start compute trace" << std::endl;
  int n = 0;
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
  for (int i=0; i < n_latent; i++) {
    VectorXd rb_trace_K (latents[i]->get_n_theta_K());
    VectorXd rb_trace_sigma (latents[i]->get_n_theta_sigma());

    // compute for K: tr(QQ^-1 dK^T diag(1/SV) K)
    for (int j=0; j < latents[i]->get_n_theta_K(); j++) {
      SparseMatrix<double> T = block_dK[i][j].transpose() * inv_SV.asDiagonal() * K;
      rb_trace_K[j] = chol_QQ.trace_num(T);
    }

    // compute for sigma: tr(Q^-1 K B_sigma.col(j)/SV K^T)
    for (int j=0; j < latents[i]->get_n_theta_sigma(); j++) {
      // build B_sigma_col_j (consider all latents)
      VectorXd BSigma_col_over_SV = VectorXd::Zero(V_sizes);
      BSigma_col_over_SV.segment(n, latents[i]->get_V_size()) = latents[i]->get_BSigma_col(j);
      BSigma_col_over_SV = BSigma_col_over_SV.cwiseProduct(inv_SV);

      SparseMatrix<double> T = K.transpose() * BSigma_col_over_SV.asDiagonal() * K;
      rb_trace_sigma[j] = chol_QQ.trace_num(T);
    }

    latents[i]->set_rb_trace(rb_trace_K, rb_trace_sigma);
    n += latents[i]->get_V_size();
  }

  // compute for theta_sigma
  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());
  for (int j=0; j < n_theta_sigma; j++) {
    SparseMatrix<double> T = A.transpose() * B_sigma.col(j).cwiseQuotient(noise_SV).asDiagonal() * A;
    rb_trace_noise_sigma[j] = chol_QQ.trace_num(T);
  }
if (debug) std::cout << "after compute trace" << std::endl;
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

// double ana_trace = chol_QQ.trace(M);
// std::cout <<"ana trace = " << ana_trace << std::endl;
// std::cout <<"num trace = " << rb_trace[j] << std::endl;
// std::cout << "-----" << std::endl;

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
  );

  return out;
}


void BlockModel::test_in_the_end() {
  sampleW_VY(); // compute QQ

  VectorXd rb_trace_K (latents[0]->get_n_theta_K());
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());

  // compute for K: tr(QQ^-1 dK^T diag(1/SV) K)
  for (int j=0; j < latents[0]->get_n_theta_K(); j++) {
    SparseMatrix<double> M = block_dK[0][j].transpose() * inv_SV.asDiagonal() * K;
    rb_trace_K[j] = chol_QQ.trace(M);
  }

  VectorXd second_term = VectorXd::Zero(latents[0]->get_n_theta_K());
  VectorXd cond_W_term = VectorXd::Zero(latents[0]->get_n_theta_K());
  for (int s=0; s < 1000; s++) {
    sampleW_VY();
    for (int j=0; j < latents[0]->get_n_theta_K(); j++) {
      second_term[j] += (block_dK[0][j].transpose() * getW()).cwiseProduct(inv_SV).dot(K * getW());
      cond_W_term[j] += (block_dK[0][j].transpose() * get_cond_W()).cwiseProduct(inv_SV).dot(K * get_cond_W());
    }

    std::cout << "diff at iter " << s << " = " << second_term/(s+1) - rb_trace_K - cond_W_term/(s+1) << std::endl;
  }
}

// marginal log likelihood  pi(Y) = pi(Y|W, V) pi(W|V) pi(V)
double BlockModel::log_likelihood() const {
  // compute log likelihood
  // pi(Y|W, noise_V) 
  VectorXd res = Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + noise_V).cwiseProduct(noise_mu);
  VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);

  // res ~ N(0, noise_SV) 
  // -0.5 * res^T noise_SV^-1 res - 0.5 logdet(noise_SV)
  double pi = atan(1)*4;
  double log_y_given_WV = -0.5 * (res.cwiseProduct(res).cwiseQuotient(noise_SV)).sum() - noise_SV.cwiseSqrt().array().log().sum() - 0.5 * n_obs * log(2*pi);
  
  // pi(block V)
  double logd_noiseV = 0;
  if (family != "normal") {
    logd_noiseV = NoiseUtil::log_density(family, noise_V, VectorXd::Ones(n_obs), B_nu, theta_nu, FALSE);
  }

  double logd_WV = 0;
  for (int i=0; i < n_latent; i++) {
    logd_WV += latents[i]->logd_W_V();
  }
  // std::cout << "res = " << res << std::endl;
  // std::cout << "logd_noiseV = " << logd_noiseV << std::endl;
  // std::cout << "logd_WV = " << logd_WV << std::endl;
  std::cout << "sigme.e" << sqrt(noise_SV(1)) << std::endl;
  std::cout << "log_y_given_WV = " << log_y_given_WV << std::endl;
  return log_y_given_WV + logd_noiseV + logd_WV;
}