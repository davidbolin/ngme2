// Implementation for block model and block_rep

#include "block.h"
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
  beta              (Rcpp::as<VectorXd>      (block_model["beta"])),
  n_obs             (Y.size()),
  n_la_params       (Rcpp::as<int>           (block_model["n_la_params"])),
  n_feff            (beta.size()),
  n_merr            (Rcpp::as<int>           (block_model["n_merr"])),
  n_repl            (Rcpp::as<int>           (block_model["n_repl"])),
  corr_measure      (Rcpp::as<bool>          (block_model["corr_measure"])),
  cor_rows          (),
  cor_cols          (),
  has_correlation   (),
  n_corr_pairs      (0),
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
  par_string        (Rcpp::as<string>     (block_model["par_string"]))
{
  // 1. Init controls
  Rcpp::List control_ngme = block_model["control_ngme"];
    const int burnin = control_ngme["burnin"];
    const double stepsize = control_ngme["stepsize"];
    bool init_sample_W = Rcpp::as<bool> (control_ngme["init_sample_W"]);
    n_gibbs     =  Rcpp::as<int>    (control_ngme["n_gibbs_samples"]);
    debug       = Rcpp::as<bool>   (control_ngme["debug"]);
    // reduce_var    =  Rcpp::as<bool>   (control_ngme["reduce_var"]);
    // reduce_power  =  Rcpp::as<double> (control_ngme["reduce_power"]);
    // threshold   =  Rcpp::as<double> (control_ngme["threshold"]);

if (debug) std::cout << "Begin Block Constructor" << std::endl;

  // 2. Init Fixed effects
  fix_flag[block_fix_beta]   = Rcpp::as<bool>        (control_ngme["fix_beta"]);
  if (beta.size() == 0) fix_flag[block_fix_beta]  = true;

  // 4. Init latent models
  Rcpp::List latents_in = block_model["models"];
  n_latent = latents_in.size(); // how many latent model
  for (int i=0; i < n_latent; ++i) {
    // construct acoording to models
    Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);
    unsigned long latent_seed = rng();
    latents.push_back(std::make_unique<Latent>(latent_in, latent_seed));
  }

if (debug) std::cout << "before set block A" << std::endl;
  /* Init A */
  int n = 0;
  for (std::vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
    setSparseBlock(&A, 0, n, (*it)->getA());
    n += (*it)->get_W_size();
  }
if (debug) std::cout << "After set block K" << std::endl;

  // 5. Init measurement noise
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

    n_nu          = (Rcpp::as<int>           (noise_in["n_nu"])),
    family = Rcpp::as<string>  (noise_in["noise_type"]);
    noise_mu = B_mu * theta_mu;
    noise_sigma = (B_sigma * theta_sigma).array().exp();
    if (family!="normal") {
      nu = Rcpp::as<double> (noise_in["nu"]);
      NoiseUtil::update_gig(family, nu, p_vec, a_vec, b_vec);
    }

  if (corr_measure) {
    cor_cols = Rcpp::as<vector<int>> (block_model["cor_cols"]);
    cor_rows = Rcpp::as<vector<int>> (block_model["cor_rows"]);
    has_correlation = Rcpp::as<vector<bool>> (block_model["has_correlation"]);
    n_corr_pairs = Rcpp::as<int> (block_model["n_corr_pairs"]);
    vector<Triplet<double>> Q_eps_triplet, dQ_eps_triplet;
    for (int i=0; i < cor_cols.size(); ++i) {
      Q_eps_triplet.push_back(Triplet<double>(cor_rows[i], cor_cols[i], cor_rows[i] == cor_cols[i]));
      if (has_correlation[cor_rows[i]]) {
        // ignore uncorrelated locations
        dQ_eps_triplet.push_back(Triplet<double>(cor_rows[i], cor_cols[i], cor_rows[i] == cor_cols[i]));
      }

      // permutation matrix for sampling correlated V
      // pmat = Rcpp::as<SparseMatrix<double>> (block_model["pmat"]);
      // pmat_inv = pmat.transpose();
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
  // if (fix_flag[block_fix_V]) var.fixV();

  // 7. Init solvers
  assemble();

  if (n_latent > 0) {
    VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    if (!corr_measure) {
      QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * A;
    }
    else{
      QQ = Q + A.transpose() * Q_eps * A;
    }
    chol_Q.analyze(Q);
    chol_QQ.analyze(QQ);
    LU_K.analyzePattern(K);
  }
if (debug) std::cout << "After init solver " << std::endl;

  // 8. optimizer related
  stepsizes = VectorXd::Constant(n_params, stepsize);
  steps_to_threshold = VectorXd::Constant(n_params, 0);
  indicate_threshold = VectorXd::Constant(n_params, 0);

if (debug) std::cout << "After init solver && before sampleW_V" << std::endl;

  if (n_latent > 0 && init_sample_W) {
    sampleW_V();
    sampleW_V();
  }
if (debug) std::cout << "Finish SampleW|V" << std::endl;
  burn_in(burnin);

if (debug) std::cout << "End Block Constructor" << std::endl;
}


void BlockModel::setW(const VectorXd& W) {
  int pos = 0;
  for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->setW(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::setPrevW(const VectorXd& W) {
  int pos = 0;
  for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_W_size();
    (*it)->setPrevW(W.segment(pos, size));
    pos += size;
  }
}

void BlockModel::setPrevV(const VectorXd& V) {
  int pos = 0;
  for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
    int size = (*it)->get_V_size();
    (*it)->setPrevV(V.segment(pos, size));
    pos += size;
  }
}

// sample W|VY
void BlockModel::sampleW_VY()
{
// if (debug) std::cout << "starting sampling W." << std::endl;
// double time = 0;
// auto timer_computeg = std::chrono::steady_clock::now();
  if (n_latent==0) return;

  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());

  // init Q and QQ
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

// time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timer_computeg).count();
// std::cout << "size of W and time of sampling is " << W.size() << " " << time << std::endl;

// if (debug) std::cout << "Finish sampling W" << std::endl;
}

// ---------------- get, set update gradient ------------------
// order is Latent, merr, feff
VectorXd BlockModel::get_parameter() {
if (debug) std::cout << "Start get_parameter"<< std::endl;
    VectorXd thetas (n_params);
    int pos = 0;
    for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
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
auto timer_computeg = std::chrono::steady_clock::now();

  VectorXd avg_gradient = VectorXd::Zero(n_params);

  // sample_uncond_noise_V()
  // sample_uncond_V();
  for (int i=0; i < n_gibbs; i++) {
// std::cout << "index of gibbs = " << i << std::endl;
    // gibbs sampling
    // order is really important here, sample W after cond V
    sample_cond_V();
    sampleW_VY();
    sample_cond_noise_V();

    VectorXd gradient = VectorXd::Zero(n_params);
    // get grad for each latent
    int pos = 0;
    for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      int theta_len = (*it)->get_n_params();
      gradient.segment(pos, theta_len) = (*it)->get_grad();
      pos += theta_len;
    }
time_compute_g += since(timer_computeg).count();

    // grad for merr
    gradient.segment(n_la_params, n_merr) = grad_theta_merr();

    // grad for fixed effects
    if (!fix_flag[block_fix_beta]) {
      gradient.segment(n_la_params + n_merr, n_feff) = (1.0/n_repl) * grad_beta();
    }
// std::cout << "here gradient = " << gradient << std::endl;
    avg_gradient += gradient;

auto timer_sampleW = std::chrono::steady_clock::now();
time_sample_w += since(timer_sampleW).count();
  }

  avg_gradient = (1.0/n_gibbs) * avg_gradient;
  gradients = avg_gradient;
  // EXAMINE the gradient to change the stepsize
  if (reduce_var) examine_gradient();

if (debug) {
  std::cout << "avg time for compute grad (ms): " << time_compute_g / n_gibbs << std::endl;
  std::cout << "avg time for sampling W(ms): " << time_sample_w / n_gibbs << std::endl;
  std::cout << "gradients = " << gradients << std::endl;
  std::cout << "Finish block gradient"<< std::endl;
}

  return gradients;
}

void BlockModel::set_parameter(const VectorXd& Theta) {
  int pos = 0;
  for (std::vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
    int theta_len = (*it)->get_n_params();
    VectorXd theta = Theta.segment(pos, theta_len);
    (*it)->set_parameter(theta);
    pos += theta_len;
  }

  // measurement noise
  set_theta_merr(Theta.segment(n_la_params, n_merr));
  if (family!="normal") NoiseUtil::update_gig(family, nu, p_vec, a_vec, b_vec);

  // fixed effects
  if (!fix_flag[block_fix_beta]) {
    beta = Theta.segment(n_la_params + n_merr, n_feff);
  }

  assemble(); //update K,dK,d2K after
  curr_iter++;
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
    W = LU_K.solve(KW);
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

  VectorXd residual = get_residual(); // + X * beta;
  VectorXd grads = X.transpose() * noise_inv_SV.asDiagonal() * residual.cwiseQuotient(noise_sigma);

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

  VectorXd residual = get_residual();
  VectorXd grad = VectorXd::Zero(n_theta_mu);
  for (int l=0; l < n_theta_mu; l++) {

      // LU_K.factorize(getK());
      // VectorXd tmp = residual + LU.solve(-h + )
      grad(l) = (noise_V - VectorXd::Ones(n_obs)).cwiseProduct(B_mu.col(l).cwiseQuotient(noise_SV)).dot(residual);
  }

  return -grad;
}

VectorXd BlockModel::grad_theta_sigma() {
  VectorXd grad = VectorXd::Zero(n_theta_sigma);
  VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);
  // grad = B_sigma.transpose() * (-0.5 * VectorXd::Ones(n_obs) + residual.array().pow(2).matrix().cwiseQuotient(noise_SV));

  VectorXd residual = get_residual();
  VectorXd vsq = (residual).array().pow(2).matrix().cwiseProduct(noise_V.cwiseInverse());
  VectorXd tmp1 = vsq.cwiseProduct(noise_sigma.array().pow(-2).matrix()) - VectorXd::Ones(n_obs);
  grad = B_sigma.transpose() * tmp1;
  // grad = - 0.5* B_sigma.transpose() * VectorXd::Ones(n_obs)
  // + B_sigma.transpose() * noise_sigma.array().pow(-2).matrix() * residual.cwiseProduct(noise_V.cwiseInverse()).dot(residual);


  // VectorXd Y_tilde_sq = (Y - A * getW() - X * beta).array().pow(2);
  // if (family == "normal") {
  //     for (int i=0; i < n_theta_sigma; i++) {
  //         grad(i) = 0.5 * B_sigma.col(i).dot(VectorXd::Ones(n_obs) - Y_tilde_sq.cwiseQuotient(noise_sigma));
  //     }
  // }

  return -grad;
}

VectorXd BlockModel::get_theta_merr() const {
  VectorXd theta_merr = VectorXd::Zero(n_merr);

  theta_merr.segment(0, n_theta_mu) = theta_mu;
  theta_merr.segment(n_theta_mu, n_theta_sigma) = theta_sigma;
  if (family != "normal")
    theta_merr(n_theta_mu + n_theta_sigma) = log(nu);

  // estimate correlation
// std::cout << " rho === " << rho << std::endl;
  if (corr_measure) theta_merr(n_merr-1) = rho2th(rho);
  // if (corr_measure) theta_merr(n_merr-1) = (rho);

  return theta_merr;
}

VectorXd BlockModel::grad_theta_merr() {
  VectorXd grad = VectorXd::Zero(n_merr);

  if (!fix_flag[block_fix_theta_mu])     grad.segment(0, n_theta_mu) = grad_theta_mu();
  if (!fix_flag[block_fix_theta_sigma])  grad.segment(n_theta_mu, n_theta_sigma) = grad_theta_sigma();
  if (!fix_flag[block_fix_nu] && family != "normal") {
    grad(n_theta_mu + n_theta_sigma) = NoiseUtil::grad_theta_nu(family, nu, noise_V, noise_prevV);
  }

  // grad of theta_rho
  if (corr_measure) {
    // Q_eps_solver.factorize(Q_eps);
    double trace = 0.5 * 2 * rho/(1-rho*rho) * n_corr_pairs;
    VectorXd res = get_residual();
    double drhs = -0.5 * (res).dot(dQ_eps * res);
    grad(n_merr-1) = trace + drhs;
    grad(n_merr-1) *= - 1.0 / (n_obs) * dtheta_th(rho);
// std::cout << "drhs = " << drhs << std::endl;
// std::cout << "trace = " << trace << std::endl;
// std::cout << "grad of rho=" << grad(n_merr-1) << std::endl;
  }

  return grad;
}

void BlockModel::set_theta_merr(const VectorXd& theta_merr) {
  theta_mu = theta_merr.segment(0, n_theta_mu);
  theta_sigma = theta_merr.segment(n_theta_mu, n_theta_sigma);
  if (family != "normal")
    nu = exp(theta_merr(n_theta_mu + n_theta_sigma));

  // update mu, sigma
  noise_mu = (B_mu * theta_mu);
  noise_sigma = (B_sigma * theta_sigma).array().exp();

  // update rho, and Q_eps
  if (corr_measure) {
    rho = th2rho(theta_merr(n_merr-1));
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
      Rcpp::Named("nu")           = nu,
      Rcpp::Named("V")            = noise_V,
      Rcpp::Named("rho")          = rho
    ),
    Rcpp::Named("beta")            = beta,
    Rcpp::Named("models")          = latents_output
  );

  return out;
}

// posterior
Rcpp::List BlockModel::sampling(int n, bool posterior) {
  std::vector<VectorXd> AWs; // blockA * blockW
  std::vector<VectorXd> Ws; // blockW
  std::vector<VectorXd> Vs; // blockV
  std::vector<VectorXd> mn_Vs; // measurement nosie V

  burn_in(5);

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
    // Vs.push_back(getV());
    // mn_Vs.push_back(var.getV());
  }

  return Rcpp::List::create(
    Rcpp::Named("AW") = AWs,
    Rcpp::Named("W") = Ws
    // ,
    // Rcpp::Named("V") = Vs,
    // Rcpp::Named("noise_V") = mn_Vs
  );
}

// fix parameter if converge
// void BlockModel::check_converge(vector<bool>& converge) {
//   int pos = 0;
//   for (std::vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
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

MatrixXd BlockModel::precond(int strategy) const {
  MatrixXd precond = MatrixXd::Zero(n_params, n_params);

  if (strategy == 0) {
    // No preconditioner
    int index_params = 0;
    for (int i=0; i < n_latent; i++) {
      int n_la = latents[i]->get_n_params();
        precond.block(index_params, index_params, n_la, n_la) = VectorXd::Constant(n_la, latents[i]->get_V_size()).asDiagonal();
        index_params += n_la;
    }

    precond.bottomRightCorner(n_merr + n_feff, n_merr + n_feff) = VectorXd::Constant(n_merr + n_feff, n_obs).asDiagonal();
  } else if (strategy == 1) {
    // Compute fast preconditioner (bdiag)
    int index_params = 0;
    for (int i=0; i < n_latent; i++) {
        MatrixXd hess = latents[i]->precond();
        precond.block(index_params, index_params, latents[i]->get_n_params(), latents[i]->get_n_params()) = latents[i]->precond();
        index_params += latents[i]->get_n_params();
    }
    // preconditioner for fixed effects and measurement error
    VectorXd v (n_merr + n_feff);
    if (family != "normal")
      v << theta_mu, theta_sigma, nu, beta;
    else
      v << theta_sigma, beta;
    precond.bottomRightCorner(n_merr + n_feff, n_merr + n_feff) = num_h_no_latent(v, 1e-5);
  } else if (strategy == 2) {
    // compute full hessian
  } else {
    throw std::invalid_argument("Invalid preconditioner strategy");
  }

// std::cout << "block precond =" << precond <<std::endl;
  return precond;
}

// precond_fast: numerical hessian for c(beta, theta_mu, theta_sigma, nu)
MatrixXd BlockModel::num_h_no_latent(const VectorXd& v, double eps) const {
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
			hessian(i, j) = (f_vij - f_v(i) - f_v(j) + original_val) / (eps * eps);
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

double BlockModel::logd_no_latent(const VectorXd& v) const {
	VectorXd theta_mu = v.segment(0, n_theta_mu);
	VectorXd theta_sigma = v.segment(n_theta_mu, n_theta_sigma);
  // n_nu can be 0 or 1
  VectorXd beta = v.segment(n_theta_mu + n_theta_sigma + n_nu, n_feff);

  // compute mu & sigma
  VectorXd noise_mu = (B_mu * theta_mu);
  VectorXd noise_sigma = (B_sigma * theta_sigma).array().exp();
// std::cout << "noise_mu=" << noise_mu << std::endl;
// std::cout << "noise_sigma=" << noise_sigma << std::endl;

  // pi(Y|W, V) (no correction for noise!)
  // residual but using prevV
	VectorXd tmp = Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + noise_prevV).cwiseProduct(noise_mu);

	double lhs = noise_sigma.array().log().sum();
	double rhs = tmp.cwiseQuotient(noise_sigma).array().square().sum();
	double logd_res = 0.5 * (lhs - rhs);

  // pi(V)
  double logd_V = 0;
  VectorXd h = VectorXd::Ones(n_obs);
  if (family != "normal") {
    double nu = v(n_theta_mu + n_theta_sigma);
    logd_V = NoiseUtil::log_density(family, noise_V, h, nu, FALSE);
  }

  return -(logd_res + logd_V);
}