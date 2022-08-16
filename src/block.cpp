#include "block.h"
#include <random>
#include <cmath>

using std::pow;

BlockModel::BlockModel(
    Rcpp::List general_in,
    Rcpp::List latents_in,
    Rcpp::List noise_in,
    Rcpp::List control_list,
    Rcpp::List debug_list
) : 
    seed          ( Rcpp::as<unsigned long> (general_in["seed"]) ),
    X             ( Rcpp::as<MatrixXd>   (general_in["X"]) ),
    Y             ( Rcpp::as<VectorXd>   (general_in["Y"]) ), 
    n_meshs       ( Rcpp::as<int>        (general_in["n_meshs"]) ),
    
    beta          ( Rcpp::as<VectorXd>   (general_in["beta"]) ),

    B_mu          ( Rcpp::as<MatrixXd>   (noise_in["B_mu"]) ),
    theta_mu      ( Rcpp::as<VectorXd>   (noise_in["theta_mu"]) ),
    n_theta_mu    ( theta_mu.size()),

    B_sigma       ( Rcpp::as<MatrixXd>   (noise_in["B_sigma"]) ),
    theta_sigma   ( Rcpp::as<VectorXd>   (noise_in["theta_sigma"]) ),
    n_theta_sigma ( theta_sigma.size()),

    n_latent      ( latents_in.size()),  // how many latent model
    
    n_obs         ( Y.size()),
    // n_params      ( Rcpp::as<int>        (general_in["n_params"]) ),
    n_la_params   ( Rcpp::as<int>        (general_in["n_la_params"]) ),
    n_feff        ( beta.size()),
    
    n_gibbs       ( Rcpp::as<int>     (control_list["gibbs_sample"]) ),
    opt_fix_effect( Rcpp::as<bool>    (control_list["opt_fix_effect"]) ),
    kill_var      ( Rcpp::as<bool>    (control_list["kill_var"]) ),
    kill_power    ( Rcpp::as<double>  (control_list["kill_power"]) ), 
    threshold     ( Rcpp::as<double>  (control_list["threshold"]) ), 
    termination   ( Rcpp::as<double>  (control_list["termination"]) ), 

    A             ( n_obs, n_meshs), 
    K             ( n_meshs, n_meshs), 
    // dK            ( n_meshs, n_meshs),
    // d2K           ( n_meshs, n_meshs),

    debug         ( Rcpp::as<bool> (debug_list["debug"]) ),
    fix_W         ( Rcpp::as<bool> (debug_list["fix_W"]) ),
    fix_merr      ( Rcpp::as<bool> (debug_list["fix_merr"])),
    residual      ( VectorXd::Zero(n_obs) )
{        
if (debug) std::cout << "Begin Block Constructor" << std::endl;  
    rng.seed(seed);
    const int burnin = control_list["burnin"];
    const double stepsize = control_list["stepsize"];

    // Init each latent model
    for (int i=0; i < n_latent; ++i) {
        Rcpp::List latent_in = Rcpp::as<Rcpp::List> (latents_in[i]);

        // construct acoording to models
        unsigned long latent_seed = rng();
        string type = latent_in["model_type"];
        if (type == "ar1") {
            latents.push_back(new AR(latent_in, latent_seed) );
        } 
        else if (type == "spde.matern") {
            latents.push_back(new Matern_ns(latent_in, latent_seed));
        } else if (type=="matern") {
            latents.push_back(new Matern(latent_in, latent_seed));
        }
    }

    // Initialize W
    // if (start_in["block.W"] != R_NilValue) {
    //     VectorXd block_W = Rcpp::as<VectorXd>   (start_in["block.W"]);
    //     // set Both W and PrevW
    //     setW(block_W); setW(block_W);
    // }

    /* Fixed effects */
    if (beta.size()==0) opt_fix_effect = false;

    /* Init variables: h, A */
    int n = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        setSparseBlock(&A,   0, n, (*it)->getA());            
        n += (*it)->getSize();
    }
    assemble();

    VectorXd inv_SV = VectorXd::Constant(n_meshs, 1).cwiseQuotient(getSV());
    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    SparseMatrix<double> QQ;
    
    /* Measurement noise */
    family = Rcpp::as<string>  (noise_in["type"]);
    noise_mu = B_mu * theta_mu;
    noise_sigma = (B_sigma * theta_sigma).array().exp();

    if (family=="normal") {
        var = new normal(n_obs);
        VectorXd theta_sigma = Rcpp::as<VectorXd>  (noise_in["theta_sigma"]);
        QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(var->getV()).asDiagonal() * A;
        n_merr = n_theta_sigma;
        // sigma_eps = theta_sigma(0);
        // QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;
    } else if (family=="nig") {
        n_merr = n_theta_sigma + n_theta_mu + 1;
        double theta_V = Rcpp::as< double >      (noise_in["theta_V"]);
        var = new ind_IG(theta_V, n_obs, rng());
    }
    
    n_params = n_la_params + n_feff + n_merr;

    chol_Q.analyze(Q);
    chol_QQ.analyze(QQ);
    LU_K.analyzePattern(K);

    // optimizer related
    stepsizes = VectorXd::Constant(n_params, stepsize);
    steps_to_threshold = VectorXd::Constant(n_params, 0);
    indicate_threshold = VectorXd::Constant(n_params, 0);

    // burn in
    sampleW_V();

    burn_in(burnin + 5);
    update_residual();
if (debug) std::cout << "init residual" << residual << std::endl;        
        
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
  for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
      int size = (*it)->getSize();
      (*it)->setW(W.segment(pos, size));
      pos += size;
  }
}

// sample W|VY 
void BlockModel::sampleW_VY()
{
// if (debug) std::cout << "starting sampling W." << std::endl;
    if (fix_W) return;

    VectorXd SV = getSV();
    VectorXd inv_SV = VectorXd::Constant(SV.size(), 1).cwiseQuotient(SV);

    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    // SparseMatrix<double> QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;
    // SparseMatrix<double> QQ = Q + A.transpose() * noise_sigma.cwiseInverse().asDiagonal() * A;
    VectorXd noise_V = var->getV();
    SparseMatrix<double> QQ = Q + A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * A;
    chol_QQ.compute(QQ);

    // VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() + 
    //     pow(sigma_eps, -2) * A.transpose() * (Y - X * beta);
    update_residual();
    VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() + 
        A.transpose() * noise_sigma.array().pow(-2).matrix().cwiseQuotient(noise_V).asDiagonal() * (residual + A * getW());
    
    VectorXd z (n_meshs); 
    z = rnorm_vec(n_meshs, 0, 1, rng());
    // sample W ~ N(QQ^-1*M, QQ^-1)
    VectorXd W = chol_QQ.rMVN(M, z);
    setW(W);

// if (debug) std::cout << "Finish sampling W" << std::endl;        
}


// ---------------- get, set update gradient ------------------
VectorXd BlockModel::get_parameter() const {
if (debug) std::cout << "Start block get parameter"<< std::endl;   
// if (debug) std::cout << "n_params = " << n_params << std::endl;   
    VectorXd thetas (n_params);
    int pos = 0;
    for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
        VectorXd theta = (*it)->get_parameter();
        thetas.segment(pos, theta.size()) = theta;
        pos += theta.size();
    }

    // fixed effects
    if (opt_fix_effect) {
        thetas.segment(n_la_params, n_feff) = beta;
    }
    
    // measurement error
    thetas.segment(n_la_params + n_feff, n_theta_sigma) = get_theta_merr();
    // thetas(n_params-1) = get_theta_sigma_eps();
if (debug) std::cout << "Theta after merr" << thetas << std::endl;   

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
        for (std::vector<Latent*>::const_iterator it = latents.begin(); it != latents.end(); it++) {
            int theta_len = (*it)->get_n_params();
            gradient.segment(pos, theta_len) = (*it)->get_grad();
            pos += theta_len;
        }
time_compute_g += since(timer_computeg).count();

        // fixed effects
        if (opt_fix_effect) {
            gradient.segment(n_la_params, n_feff) = grad_beta();
        }
if (debug) std::cout << "Finish fixed effects"<< std::endl;   
        
        // gradient.segment(n_la_params + n_feff, n_theta_sigma) = grad_theta_merr();
        gradient.segment(n_la_params + n_feff, n_merr) = grad_theta_merr();
if (debug) std::cout << "Finish measurement noise"<< std::endl;   

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

if (debug) std::cout << "Finish block gradient"<< std::endl;   
    return gradients;
}


void BlockModel::set_parameter(const VectorXd& Theta) {
    int pos = 0;
    for (std::vector<Latent*>::iterator it = latents.begin(); it != latents.end(); it++) {
        int theta_len = (*it)->get_n_params();
        VectorXd theta = Theta.segment(pos, theta_len);
        (*it)->set_parameter(theta);
        pos += theta_len;
    }

    // fixed effects
    if (opt_fix_effect) {
        beta = Theta.segment(n_la_params, n_feff);
    }

    // measurement noise
    set_theta_merr(Theta.segment(n_la_params + n_feff, n_theta_sigma));
    // set_theta_sgima_eps(Theta(n_params-1));
    
    assemble(); //update K,dK,d2K after
    update_residual();
// std::cout << "Theta=" << Theta <<std::endl;
}

// sample W|V
void BlockModel::sampleW_V()
{
  // sample KW ~ N(mu*(V-h), diag(V))
  VectorXd SV = getSV();
  Eigen::VectorXd KW (n_meshs);
  for (int i=0; i < n_meshs; i++) {
    KW[i] = R::rnorm(0, sqrt(SV[i]));
  }
  KW = getMean() + KW;

  LU_K.factorize(K);
  VectorXd W = LU_K.solve(KW);

  setW(W);
}


// --------- Measurement Noise ---------------
VectorXd BlockModel::grad_beta() {
    VectorXd noise_V = var->getV();
    VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());
    
    update_residual();
    VectorXd grads = X.transpose() * noise_inv_SV.asDiagonal() * residual;
    MatrixXd hess = X.transpose() * noise_inv_SV.asDiagonal() * X;
    grads = hess.ldlt().solve(grads);

// std::cout << "grads of beta=" << -grads << std::endl;        
    return -grads;
}

VectorXd BlockModel::grad_theta_mu() {
    VectorXd noise_V = var->getV();
    VectorXd noise_inv_SV = noise_V.cwiseProduct(noise_sigma.array().pow(-2).matrix());
    MatrixXd noise_X = (-VectorXd::Ones(n_obs) + noise_V).asDiagonal() * B_mu;
    update_residual();
    VectorXd grads = noise_X.transpose() * noise_inv_SV.asDiagonal() * residual;
    
    MatrixXd hess = noise_X.transpose() * noise_inv_SV.asDiagonal() * noise_X;
    grads = hess.ldlt().solve(grads);

    return -grads;
}

VectorXd BlockModel::grad_theta_sigma() {
    VectorXd grad = VectorXd::Zero(n_theta_sigma);
// std::cout << "noise_sigma =" << noise_sigma << std::endl;        
// std::cout << "residual =" << residual << std::endl;        
// std::cout << "th_sigma =" << residual << std::endl;        

std::cout << "print W =" << getW() << std::endl; 
    VectorXd noise_V = var->getV();
    VectorXd noise_SV = noise_sigma.array().pow(2).matrix().cwiseProduct(noise_V);
    // grad = B_sigma.transpose() * (-0.5 * VectorXd::Ones(n_obs) + residual.array().pow(2).matrix().cwiseQuotient(noise_SV));
    
    update_residual();
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
std::cout << "grads of th_sigma =" << grad << std::endl;        
    return grad;
}

void BlockModel::update_residual() {
    residual = Y - A * getW() - X * beta - (-VectorXd::Ones(n_obs) + var->getV()).cwiseProduct(noise_mu);
}
