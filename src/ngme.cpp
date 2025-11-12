// implement the Ngme class and rand effect class
#include "ngme.h"

#ifdef _OPENMP
  #include<omp.h>
  #pragma omp declare reduction(vec_plus : Eigen::VectorXd : omp_out += omp_in) initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))
  #pragma omp declare reduction(mat_plus : Eigen::MatrixXd : omp_out += omp_in) initializer(omp_priv = Eigen::MatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))
#endif

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List& R_ngme, unsigned long seed, int sampling_strategy, int num_threads_repl, double start_sd) :
  n_params          (Rcpp::as<int> (R_ngme["n_params"])),
  sampling_strategy (sampling_strategy),
  num_threads_repl  (num_threads_repl),
  gen               (seed),
  debug             (false)
{
  Rcpp::List control_ngme = R_ngme["control_ngme"];
    debug = Rcpp::as<bool> (control_ngme["debug"]);
    stepsize = VectorXd::Constant(n_params, Rcpp::as<double> (control_ngme["stepsize"]));

  Rcpp::List list_replicates = Rcpp::as<Rcpp::List> (R_ngme["replicates"]);
  n_repl = list_replicates.size();

  // build for each replicates
  for (int i=0; i < n_repl; i++) {
    Rcpp::List ngme_repl = Rcpp::as<Rcpp::List> (list_replicates[i]);
    ngme_repls.push_back(std::make_shared<BlockModel>(ngme_repl, seed + i));
    num_each_repl.push_back(ngme_repls[i]->get_n_obs());
  }
  sum_num_each_repl = std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);

  // Init the random number generator
  weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(), num_each_repl.end());

  VectorXd p = ngme_repls[0]->get_parameter();
  std::normal_distribution<double> distribution(0.0, start_sd);
  for (int i = 0; i < p.size(); ++i) {
      p[i] += distribution(gen);
  }
  ngme_repls[0]->set_parameter(p);
}

void Ngme::compute(bool with_precond, double eps) {
  // Aggregate over replicates
  last_grad_ = VectorXd::Zero(n_params);
  // If with_precond=true, also compute and cache preconditioner for current strategy
  if (with_precond) {
    last_precond_ = MatrixXd::Zero(n_params, n_params);
    precond_valid_ = false;
  }
  if (sampling_strategy == Strategy::all) {
    for (int i=0; i < n_repl; i++) {
      if (with_precond) ngme_repls[i]->set_precond_strategy(curr_precond_strategy_);
      ngme_repls[i]->compute_grad_and_hessian(with_precond, eps);
      last_grad_ += ngme_repls[i]->get_grad_with_gibbs_samples();
      if (with_precond) last_precond_ += ngme_repls[i]->get_preconditioner(curr_precond_strategy_, eps);
    }
  } else { // ws
    int idx = weighted_sampler(gen);
    if (with_precond) ngme_repls[idx]->set_precond_strategy(curr_precond_strategy_);
    ngme_repls[idx]->compute_grad_and_hessian(with_precond, eps);
    last_grad_ = ngme_repls[idx]->get_grad_with_gibbs_samples();
    if (with_precond) last_precond_ = ngme_repls[idx]->get_preconditioner(curr_precond_strategy_, eps);
  }
  grad_valid_ = true;
  if (with_precond) {
    last_precond_eps_ = eps;
    last_precond_strategy_ = curr_precond_strategy_;
    precond_valid_ = true;
  }
}

MatrixXd Ngme::precond(double eps) {
  // Pure getter: return cached preconditioner if computed for current strategy
  if (!precond_valid_ || std::abs(eps - last_precond_eps_) > 1e-12) {
    // fallback: identity-like
    return VectorXd::Constant(n_params, 1e-5).asDiagonal();
  }
  return last_precond_;
}

VectorXd Ngme::grad() {
  if (!grad_valid_) {
    // fallback: compute gradients only
    compute(false, 1e-5);
  }
  return last_grad_;
}


void Ngme::burn_in(int iterations) {
  #pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i=0; i < n_repl; i++) {
    ngme_repls[i]->burn_in(iterations);
  }
}

VectorXd Ngme::get_parameter() {
  VectorXd p = ngme_repls[0]->get_parameter();
  return p;
}

void Ngme::set_parameter(const VectorXd& p) {
  // set the same parameter for all latent
// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  #pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i=0; i < n_repl; i++) {
    ngme_repls[i]->set_parameter(p);
  }

  // set the different parameter for each random effect
  if (debug) std::cout << "set_parameter() in ngme class" << std::endl;
  // Invalidate caches
  grad_valid_ = false;
  precond_valid_ = false;
}
