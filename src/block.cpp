#include "block.h"
#include <random>
#include <cmath>

using std::pow;

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
  n_params          (Rcpp::as<int>           (block_model["n_params"])),
  n_feff            (beta.size()),
  n_merr            (Rcpp::as<int>           (block_model["n_merr"])),

  debug             (Rcpp::as<bool>          (block_model["debug"])),
  A                 (n_obs, W_sizes),
  K                 (V_sizes, W_sizes),
  var               (Var(Rcpp::as<Rcpp::List> (block_model["noise"]), rng())),

  curr_iter         (0),
  beta_traj         (beta.size()),
  // dK            (V_sizes, W_sizes)
  // d2K           (V_sizes, W_sizes)
  par_string        (Rcpp::as<string>     (block_model["par_string"]))
{
  // 1. Init controls
  Rcpp::List control_in = block_model["control"];
    const int burnin = control_in["burnin"];
    const double stepsize = control_in["stepsize"];
    n_gibbs     =  Rcpp::as<int>    (control_in["gibbs_sample"]);
    opt_beta    =  Rcpp::as<bool>   (control_in["opt_beta"]);
    kill_var    =  Rcpp::as<bool>   (control_in["kill_var"]);
    kill_power  =  Rcpp::as<double> (control_in["kill_power"]);
    threshold   =  Rcpp::as<double> (control_in["threshold"]);
    termination =  Rcpp::as<double> (control_in["termination"]);

if (debug) std::cout << "Begin Block Constructor" << std::endl;

  // 2. Init Fixed effects
  fix_flag[block_fix_beta]   = Rcpp::as<bool>        (control_in["fix_beta"]);
  if (beta.size() == 0) opt_beta = false;

  // 3. Init latent models
  Rcpp::List latents_in = block_model["latents"];
  n_latent = latents_in.size(); // how many latent model
  for (int i=0; i < n_latent; ++i) {
    Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);
    int n_theta_K = Rcpp::as<int> (latent_in["n_theta_K"]);

    // construct acoording to models
    unsigned long latent_seed = rng();
    string model_type = latent_in["model"];
    if (model_type == "ar1") {
      latents.push_back(std::make_unique<AR>(latent_in, latent_seed));
    }
    else if (model_type == "rw1") {
      latents.push_back(std::make_unique<AR>(latent_in, latent_seed));
    }
    else if (model_type == "matern" && n_theta_K > 1) {
      latents.push_back(std::make_unique<Matern_ns>(latent_in, latent_seed));
    } else if (model_type=="matern" && n_theta_K == 1) {
      latents.push_back(std::make_unique<Matern>(latent_in, latent_seed));
    } else {
      std::cout << "Unknown model." << std::endl;
    }
  }

  /* Init variables: h, A */
  int n = 0;
  for (std::vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
    setSparseBlock(&A,   0, n, (*it)->getA());
    n += (*it)->get_W_size();
  }
  assemble();

if (debug) std::cout << "After block assemble" << std::endl;

  // 4. Init measurement noise
  Rcpp::List noise_in   = block_model["noise"];

  B_mu          = (Rcpp::as<MatrixXd>      (noise_in["B_mu"])),
  theta_mu      = (Rcpp::as<VectorXd>      (noise_in["theta_mu"])),
  n_theta_mu    = (theta_mu.size()),

  B_sigma       = (Rcpp::as<MatrixXd>      (noise_in["B_sigma"])),
  theta_sigma   = (Rcpp::as<VectorXd>      (noise_in["theta_sigma"])),
  n_theta_sigma = (theta_sigma.size()),

  fix_flag[block_fix_theta_mu]        = Rcpp::as<bool> (noise_in["fix_theta_mu"]);
  fix_flag[block_fix_theta_sigma]     = Rcpp::as<bool> (noise_in["fix_theta_sigma"]);

  family = Rcpp::as<string>  (noise_in["noise_type"]);
  noise_mu = B_mu * theta_mu;
  noise_sigma = (B_sigma * theta_sigma).array().exp();

if (debug) std::cout << "After block construct noise" << std::endl;

  // 5. Fix V and init V
  // if (fix_flag[block_fix_V]) var.fixV();

  // 6. Init solvers
  if(n_latent > 0){
    VectorXd inv_SV = VectorXd::Constant(V_sizes, 1).cwiseQuotient(getSV());
    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    SparseMatrix<double> QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(var.getV()).asDiagonal() * A;
    chol_Q.analyze(Q);
    chol_QQ.analyze(QQ);
    LU_K.analyzePattern(K);
  }

  // 7. optimizer related
  stepsizes = VectorXd::Constant(n_params, stepsize);
  steps_to_threshold = VectorXd::Constant(n_params, 0);
  indicate_threshold = VectorXd::Constant(n_params, 0);

  // 8. Do the burn-in
  if(n_latent > 0) {
    sampleW_V();
    sampleW_V();
  }

  // record
  theta_mu_traj.resize(n_theta_mu);
  theta_sigma_traj.resize(n_theta_sigma);
  record_traj();

if (debug) std::cout << "End Block Constructor" << std::endl;
}


// ---- helper function for sampleW ----
Eigen::VectorXd rnorm_vec(int n, double mu, double sigma, unsigned long seed=0)
{
  std::mt19937 norm_rng(seed);
  std::normal_distribution<double> rnorm {0,1};
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; i++)
  {
    // out[i] = R::rnorm(mu, sigma);
    out[i] = rnorm(norm_rng) * sigma + mu;
  }
  return (out);
}

// ---- other functions ------
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
  if (n_latent==0) return;

  VectorXd SV = getSV();
  VectorXd inv_SV = VectorXd::Constant(SV.size(), 1).cwiseQuotient(SV);
  // VectorXd V = getV();
  // VectorXd inv_V = VectorXd::Constant(V.size(), 1).cwiseQuotient(V);

  SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
  // SparseMatrix<double> QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;
  // SparseMatrix<double> QQ = Q + A.transpose() * noise_sigma.cwiseInverse().asDiagonal() * A;
  VectorXd noise_V = var.getV();
  SparseMatrix<double> QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * A;
  chol_QQ.compute(QQ);

  // VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() +
  //     pow(sigma_eps, -2) * A.transpose() * (Y - X * beta);
  VectorXd residual = get_residual();
  VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() +
  // VectorXd M = K.transpose() * inv_V.asDiagonal() * getMean() +
      A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * (residual + A * getW());

  VectorXd z (W_sizes);
  z = rnorm_vec(W_sizes, 0, 1, rng());
  // sample W ~ N(QQ^-1*M, QQ^-1)
  VectorXd W = chol_QQ.rMVN(M, z);
  setW(W);

// if (debug) std::cout << "Finish sampling W" << std::endl;
}

// ---------------- get, set update gradient ------------------
VectorXd BlockModel::get_parameter() const {
if (debug) std::cout << "Start block get parameter"<< std::endl;
    VectorXd thetas (n_params);
    int pos = 0;
    for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      VectorXd theta = (*it)->get_parameter();
      thetas.segment(pos, theta.size()) = theta;
      pos += theta.size();
    }

    if (opt_beta) {
      thetas.segment(n_la_params, n_feff) = beta;
    }

    thetas.segment(n_la_params + n_feff, n_merr) = get_theta_merr();

if (debug) std::cout << "Finish block get parameter"<< std::endl;
    return thetas;
}


VectorXd BlockModel::grad() {
if (debug) std::cout << "Start block gradient"<< std::endl;
long long time_compute_g = 0;
long long time_sample_w = 0;

  VectorXd avg_gradient = VectorXd::Zero(n_params);
  for (int i=0; i < n_gibbs; i++) {
    // stack grad
    VectorXd gradient = VectorXd::Zero(n_params);

auto timer_computeg = std::chrono::steady_clock::now();
    // get grad for each latent
    int pos = 0;
    for (std::vector<std::unique_ptr<Latent>>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      int theta_len = (*it)->get_n_params();
      gradient.segment(pos, theta_len) = (*it)->get_grad();
      pos += theta_len;
    }
time_compute_g += since(timer_computeg).count();

    // fixed effects
    if (opt_beta) {
      gradient.segment(n_la_params, n_feff) = grad_beta();
    }

    // gradient.segment(n_la_params + n_feff, n_theta_sigma) = grad_theta_merr();
    gradient.segment(n_la_params + n_feff, n_merr) = grad_theta_merr();

    avg_gradient += gradient;

    // gibbs sampling
    sampleV_WY();
auto timer_sampleW = std::chrono::steady_clock::now();
    sampleW_VY();
time_sample_w += since(timer_sampleW).count();
    sample_cond_block_V();
  }

if (debug) {
std::cout << "avg time for compute grad (ms): " << time_compute_g / n_gibbs << std::endl;
std::cout << "avg time for sampling W(ms): " << time_sample_w / n_gibbs << std::endl;
}

  avg_gradient = (1.0/n_gibbs) * avg_gradient;
  gradients = avg_gradient;
  // EXAMINE the gradient to change the stepsize
  if (kill_var) examine_gradient();

if (debug) std::cout << "gradients = " << gradients << std::endl;
if (debug) std::cout << "Finish block gradient"<< std::endl;
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

  // fixed effects
  if (opt_beta) {
    beta = Theta.segment(n_la_params, n_feff);
  }

    // measurement noise
  set_theta_merr(Theta.segment(n_la_params + n_feff, n_merr));
  record_traj();

  assemble(); //update K,dK,d2K after
  curr_iter++;
}

// sample W|V
void BlockModel::sampleW_V()
{
  if(n_latent==0) return;
  std::normal_distribution<double> rnorm {0,1};

  // sample KW ~ N(mu*(V-h), diag(V))
  VectorXd SV = getSV();
  Eigen::VectorXd KW (V_sizes);
  for (int i=0; i < V_sizes; i++) {
    // KW[i] = R::rnorm(0, sqrt(SV[i]));
    KW[i] = rnorm(rng) * sqrt(SV[i]);
  }
  KW = getMean() + KW;

  VectorXd W (W_sizes);
  if (V_sizes == W_sizes) {
    LU_K.factorize(K);
    W = LU_K.solve(KW);
  } else {
    SparseMatrix<double> Q = K.transpose() * K;
    chol_Q.compute(Q);
    W = chol_Q.solve(K.transpose() * KW);
  }

  setW(W);
}

// --------- Fiexed effects and Measurement Error ---------------
VectorXd BlockModel::grad_beta() {
  VectorXd noise_V = var.getV();
  VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());

  VectorXd residual = get_residual();
  VectorXd grads = X.transpose() * noise_inv_SV.asDiagonal() * residual;
  MatrixXd hess = X.transpose() * noise_inv_SV.asDiagonal() * X;
  grads = hess.ldlt().solve(grads);
//   Rcpp::Rcout << "(beta) grads = " << grads << "\n";
// std::cout << "grads of beta=" << -grads << std::endl;
    return -grads;
}

VectorXd BlockModel::grad_theta_mu() {
  // VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());
  // MatrixXd noise_X = (-VectorXd::Ones(n_obs) + noise_V).asDiagonal() * B_mu;

  // VectorXd grad = noise_X.transpose() * noise_inv_SV.asDiagonal() * residual;
  // MatrixXd hess = noise_X.transpose() * noise_inv_SV.asDiagonal() * noise_X;
  // grad = hess.ldlt().solve(grad);

  VectorXd noise_V = var.getV();
  VectorXd noise_SV = noise_V.cwiseProduct(noise_sigma.array().pow(2).matrix());

  VectorXd residual = get_residual();
  VectorXd grad (n_theta_mu);
  for (int l=0; l < n_theta_mu; l++) {

      // LU_K.factorize(getK());
      // VectorXd tmp = residual + LU.solve(-h + )
      grad(l) = (noise_V - VectorXd::Ones(n_obs)).cwiseProduct(B_mu.col(l).cwiseQuotient(noise_SV)).dot(residual);
  }
  grad = - 1.0 / n_obs * grad;
  return grad;
}

VectorXd BlockModel::grad_theta_sigma() {
  VectorXd grad = VectorXd::Zero(n_theta_sigma);
  VectorXd noise_V = var.getV();
  VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);
  // grad = B_sigma.transpose() * (-0.5 * VectorXd::Ones(n_obs) + residual.array().pow(2).matrix().cwiseQuotient(noise_SV));

  VectorXd residual = get_residual();
  VectorXd vsq = (residual).array().pow(2).matrix().cwiseProduct(noise_V.cwiseInverse());
  VectorXd tmp1 = vsq.cwiseProduct(noise_sigma.array().pow(-2).matrix()) - VectorXd::Constant(n_obs, 1);
  grad = B_sigma.transpose() * tmp1;

  // grad = - 0.5* B_sigma.transpose() * VectorXd::Ones(n_obs)
  // + B_sigma.transpose() * noise_sigma.array().pow(-2).matrix() * residual.cwiseProduct(noise_V.cwiseInverse()).dot(residual);


  // VectorXd Y_tilde_sq = (Y - A * getW() - X * beta).array().pow(2);
  // if (family == "normal") {
  //     for (int i=0; i < n_theta_sigma; i++) {
  //         grad(i) = 0.5 * B_sigma.col(i).dot(VectorXd::Ones(n_obs) - Y_tilde_sq.cwiseQuotient(noise_sigma));
  //     }
  // }
  grad = - grad / n_obs;
  return grad;
}

VectorXd BlockModel::get_theta_merr() const {
  VectorXd theta_merr = VectorXd::Zero(n_merr);

  if (family=="normal") {
      theta_merr = theta_sigma;
  } else {
      theta_merr.segment(0, n_theta_mu) = theta_mu;
      theta_merr.segment(n_theta_mu, n_theta_sigma) = theta_sigma;
      theta_merr(n_theta_mu + n_theta_sigma) =  var.get_unbound_theta_V();
  }

  return theta_merr;
}

VectorXd BlockModel::grad_theta_merr() {
  VectorXd grad = VectorXd::Zero(n_merr);

  if (family=="normal") {
    if (!fix_flag[block_fix_theta_sigma])  grad = grad_theta_sigma();
  } else {
    if (!fix_flag[block_fix_theta_mu])     grad.segment(0, n_theta_mu) = grad_theta_mu();
    if (!fix_flag[block_fix_theta_sigma])  grad.segment(n_theta_mu, n_theta_sigma) = grad_theta_sigma();
    grad(n_theta_mu + n_theta_sigma) = var.grad_theta_var();
  }

  return grad;
}

void BlockModel::set_theta_merr(const VectorXd& theta_merr) {
  if (family=="normal") {
    theta_sigma = theta_merr;
  } else {
    theta_mu = theta_merr.segment(0, n_theta_mu);
    theta_sigma = theta_merr.segment(n_theta_mu, n_theta_sigma);
    var.set_theta_var(theta_merr(n_theta_mu + n_theta_sigma));
  }

  // update mu, sigma
  noise_mu = (B_mu * theta_mu);
  noise_sigma = (B_sigma * theta_sigma).array().exp();
}

// generate output to R
Rcpp::List BlockModel::output() const {
  Rcpp::List latents_output;
  for (int i=0; i < n_latent; i++) {
    latents_output.push_back((*latents[i]).output());
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("noise")     = Rcpp::List::create(
      Rcpp::Named("noise_type")   = family,
      Rcpp::Named("theta_mu")     = theta_mu,
      Rcpp::Named("theta_sigma")  = theta_sigma,
      Rcpp::Named("theta_V")      = var.get_theta_V(),
      Rcpp::Named("V")            = var.getV()
    ),
    Rcpp::Named("beta")             = beta,
    Rcpp::Named("latents")          = latents_output
  );

  out.attr("trajectory") = Rcpp::List::create(
    Rcpp::Named("beta")        = beta_traj,
    Rcpp::Named("theta_mu")    = theta_mu_traj,
    Rcpp::Named("theta_sigma") = theta_sigma_traj,
    Rcpp::Named("theta_V")     = theta_V_traj
  );
  return out;
}

// posterior
Rcpp::List BlockModel::sampling(int iterations, bool posterior) {
  std::vector<VectorXd> Ws;
  std::vector<VectorXd> Vs;
  std::vector<VectorXd> Block_Vs;

  burn_in(5);

  for (int i=0; i < iterations; i++) {
    if (posterior) {
      sampleV_WY();
      sampleW_VY();
      sample_cond_block_V();
    } else {
      sample_V();
      sampleW_V();
      var.sample_V();
      // construct the Y in R
    }

    Ws.push_back(getW());
    Vs.push_back(getV());
    Block_Vs.push_back(var.getV());
  }

  return Rcpp::List::create(
    Rcpp::Named("Ws") = Ws,
    Rcpp::Named("Vs") = Vs,
    Rcpp::Named("Block_Vs") = Block_Vs
  );
}

// provide stepsize
inline void BlockModel::examine_gradient() {

    // examine if the gradient under the threshold
    // for (int i=0; i < n_params; i++) {
    //     if (abs(gradients(i)) < threshold) {
    //         indicate_threshold(i) = 1;
    //     }      // mark if under threshold
    //     if (!indicate_threshold(i)) steps_to_threshold(i) = counting; // counting up
    // }

    // counting += 1;
    // stepsizes = (VectorXd::Constant(n_params, counting) - steps_to_threshold).cwiseInverse().array().pow(kill_power);

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
