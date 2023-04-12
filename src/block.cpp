// Implementation for block model and block_rep

#include "block.h"
#include <random>
#include <cmath>

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
  n_params          (Rcpp::as<int>           (block_model["n_params"])),
  n_la_params       (Rcpp::as<int>           (block_model["n_la_params"])),
  n_feff            (beta.size()),
  n_merr            (Rcpp::as<int>           (block_model["n_merr"])),
  // number of total replicates(blocks)
  n_repl            (Rcpp::as<int>           (block_model["n_repl"])),

  debug             (false),
  A                 (n_obs, W_sizes),
  K                 (V_sizes, W_sizes),
  Q                 (W_sizes, W_sizes),
  QQ                (W_sizes, W_sizes),
  var               (Var(Rcpp::as<Rcpp::List> (block_model["noise"]), rng())),

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
  Rcpp::List latents_in = block_model["latents"];
  n_latent = latents_in.size(); // how many latent model
  for (int i=0; i < n_latent; ++i) {
    // construct acoording to models
    Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);
    unsigned long latent_seed = rng();
    latents.push_back(LatentFactory::create(latent_in, latent_seed));
  }

if (debug) std::cout << "before set block A" << std::endl;
  /* Init A */
  int n = 0;
  for (std::vector<std::unique_ptr<Latent>>::iterator it = latents.begin(); it != latents.end(); it++) {
    setSparseBlock(&A, 0, n, (*it)->getA());
    n += (*it)->get_W_size();
  }

// if (debug) std::cout << "After set block K" << std::endl;

  // 5. Init measurement noise
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

  // 6. Fix V and init V
  // if (fix_flag[block_fix_V]) var.fixV();

  // 7. Init solvers
  assemble();
  if(n_latent > 0) {
    VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(getSV());
    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    SparseMatrix<double> QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(var.getV()).asDiagonal() * A;
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

// sample (U,W)|VY
void BlockModel::sampleW_VY()
{
// if (debug) std::cout << "starting sampling W." << std::endl;
  if (n_latent==0) return;

  VectorXd SV = getSV();
  VectorXd inv_SV = VectorXd::Ones(V_sizes).cwiseQuotient(SV);
  VectorXd noise_V = var.getV();

  // init Q and QQ
  Q = K.transpose() * inv_SV.asDiagonal() * K;
// std::cout << "Q = " << Q << std::endl;
  QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * A;
// std::cout << "QQ = " << QQ << std::endl;
  VectorXd residual = get_residual();
  VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean();
  // concat Solve(Q, b) for effects and process

  M += A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * (residual + A * getW());

// std::cout << "M = " << M << std::endl;
  VectorXd z (W_sizes);
  z = rnorm_vec(W_sizes, 0, 1, rng());

  // sample UW ~ N(QQ^-1*M, QQ^-1)
  chol_QQ.compute(QQ);
  VectorXd W = chol_QQ.rMVN(M, z);
  setW(W);

if (debug) std::cout << "Finish sampling W" << std::endl;
}

// ---------------- get, set update gradient ------------------
VectorXd BlockModel::get_parameter() const {
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
VectorXd BlockModel::precond_grad() {
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

    // gradient.segment(n_la_params + n_feff, n_theta_sigma) = grad_theta_merr();
    gradient.segment(n_la_params, n_merr) = grad_theta_merr();
    // fixed effects
    if (!fix_flag[block_fix_beta]) {
      gradient.segment(n_la_params + n_merr, n_feff) = (1.0/n_repl) * grad_beta();
    }

    avg_gradient += gradient;

    // gibbs sampling
    sampleV_WY();
auto timer_sampleW = std::chrono::steady_clock::now();
    sampleW_VY();
time_sample_w += since(timer_sampleW).count();
    sample_cond_block_V();
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
// std::cout << "WWWW = " << W << std::endl;
  setW(W);
}

// --------- Fiexed effects and Measurement Error ---------------
VectorXd BlockModel::grad_beta() {
  VectorXd noise_V = var.getV();
  VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());

  VectorXd residual = get_residual(); // + X * beta;
  VectorXd grads = X.transpose() * noise_inv_SV.asDiagonal() * residual.cwiseQuotient(noise_sigma);
  MatrixXd hess = X.transpose() * noise_inv_SV.asDiagonal() * X;
  // grads = grads / A.rows();
  grads = hess.ldlt().solve(grads);
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
  VectorXd grad = VectorXd::Zero(n_theta_mu);
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
      theta_merr(n_theta_mu + n_theta_sigma) =  var.get_log_nu();
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
    grad(n_theta_mu + n_theta_sigma) = var.grad_log_nu();
  }

  return grad;
}

void BlockModel::set_theta_merr(const VectorXd& theta_merr) {
  if (family=="normal") {
    theta_sigma = theta_merr;
  } else {
    theta_mu = theta_merr.segment(0, n_theta_mu);
    theta_sigma = theta_merr.segment(n_theta_mu, n_theta_sigma);
    var.set_log_nu(theta_merr(n_theta_mu + n_theta_sigma));
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
      Rcpp::Named("nu")      = var.get_nu(),
      Rcpp::Named("V")            = var.getV()
    ),
    Rcpp::Named("beta")             = beta,
    Rcpp::Named("latents")          = latents_output
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
      sampleV_WY();
      sampleW_VY();
      sample_cond_block_V();
    } else {
      sample_V();
      sampleW_V();
      var.sample_V();
      // construct the Y in R
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
inline void BlockModel::examine_gradient() {

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