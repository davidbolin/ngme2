#include "block.h"
#include <cmath>

using std::pow;

// ---- helper function for sampleW ----
Eigen::VectorXd rnorm_vec(int n, double mu, double sigma)
{
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; i++)
  {
    out[i] = R::rnorm(mu, sigma);
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
if (debug) std::cout << "starting sampling W." << std::endl;

  if (fix_W) return;

  if (family=="normal") {
    VectorXd SV = getSV();
    VectorXd inv_SV = VectorXd::Constant(SV.size(), 1).cwiseQuotient(SV);

    SparseMatrix<double> Q = K.transpose() * inv_SV.asDiagonal() * K;
    SparseMatrix<double> QQ = Q + pow(sigma_eps, -2) * A.transpose() * A;
    chol_QQ.compute(QQ);

    VectorXd M = K.transpose() * inv_SV.asDiagonal() * getMean() + 
        pow(sigma_eps, -2) * A.transpose() * (Y - X * beta);

    VectorXd z (n_meshs); 
    z = rnorm_vec(n_meshs, 0, 1);
    
    // sample W ~ N(QQ^-1*M, QQ^-1)
    VectorXd W = chol_QQ.rMVN(M, z);
    setW(W);
  } 
  else if (family=="nig") {
    // update W here
    
  }

if (debug) std::cout << "Finish sampling W" << std::endl;        
}



// ---------------- get, set update gradient ------------------
VectorXd BlockModel::get_parameter() const {
if (debug) std::cout << "Start block get parameter"<< std::endl;   
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
    thetas.segment(n_la_params + n_feff, n_merr) = get_theta_merr();
    // thetas(n_params-1) = get_theta_sigma_eps();

if (debug) std::cout << "Finish block get parameter"<< std::endl;   
    return thetas;
}


VectorXd BlockModel::grad() {
if (debug) std::cout << "Start block gradient"<< std::endl;   
    VectorXd avg_gradient = VectorXd::Zero(n_params);
    
long long time_compute_g = 0;
long long time_sample_w = 0;

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
        
        // measurement error 
        gradient.segment(n_la_params + n_feff, n_merr) = grad_theta_merr();

        avg_gradient += gradient;

        // gibbs sampling
        sampleV_WY(); 
auto timer_sampleW = std::chrono::steady_clock::now();
        sampleW_VY();
time_sample_w += since(timer_sampleW).count();
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
    set_theta_merr(Theta.segment(n_la_params + n_feff, n_merr));
    // set_theta_sgima_eps(Theta(n_params-1));
    
    assemble(); //update K,dK,d2K after
// std::cout << "Theta=" << Theta <<std::endl;
}




// sample W|V
// void BlockModel::sampleW_V()
// {
//   // sample KW ~ N(mu*(V-h), diag(V))
//   VectorXd SV = getSV();
//   Eigen::VectorXd KW (n_meshs);
//   for (int i=0; i < n_meshs; i++) {
//     KW[i] = R::rnorm(0, sqrt(SV[i]));
//   }
//   KW = getMean() + KW;

//   LU_K.factorize(K);
//   VectorXd W = LU_K.solve(KW);

//   setW(W);
// }


