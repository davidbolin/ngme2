// implement the Ngme class and rand effect calss
#include "ngme.h"

#ifdef _OPENMP
  #include<omp.h>
  #pragma omp declare reduction(vec_plus : Eigen::VectorXd : omp_out += omp_in) initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))
  #pragma omp declare reduction(mat_plus : Eigen::MatrixXd : omp_out += omp_in) initializer(omp_priv = Eigen::MatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))
#endif

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List& R_ngme, unsigned long seed, int sampling_strategy, int num_threads_repl) :
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

  // build for each replciates
  for (int i=0; i < n_repl; i++) {
    Rcpp::List ngme_repl = Rcpp::as<Rcpp::List> (list_replicates[i]);
    ngme_repls.push_back(std::make_shared<BlockModel>(ngme_repl, seed));
    num_each_repl.push_back(ngme_repls[i]->get_n_obs());
  }
  sum_num_each_repl = std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);

  // Init the random number generator
  weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(), num_each_repl.end());
}

MatrixXd Ngme::precond(int strategy, double eps) {
  MatrixXd precond = MatrixXd::Zero(n_params, n_params);

  if (sampling_strategy == Strategy::all) {
    #pragma omp parallel for schedule(static) reduction(mat_plus:precond) num_threads(num_threads_repl)
    for (int i=0; i < n_repl; i++) {
      precond += ngme_repls[i]->precond(strategy, eps) / n_repl;
    }
  } else if (sampling_strategy == Strategy::ws) {
    // weighted sampling (WS) for each replicate
    int idx = weighted_sampler(gen);
    precond = (num_each_repl[idx] / sum_num_each_repl) * ngme_repls[idx]->precond(strategy, eps);
  }

// std::cout << "precond in ngme class = \n" << precond << std::endl;
  return precond;
}

VectorXd Ngme::grad() {
  VectorXd g = VectorXd::Zero(n_params);
  // weighted averge over all replicates
  if (sampling_strategy == Strategy::all) {
    #pragma omp parallel for schedule(static) reduction(vec_plus:g) num_threads(num_threads_repl)
    for (int i=0; i < n_repl; i++) {
      g +=  (num_each_repl[i] / sum_num_each_repl) * ngme_repls[i]->grad() / n_repl;
    }
  } else if (sampling_strategy == Strategy::ws) {
    // weighted sampling (WS) for each replicate
    int idx = weighted_sampler(gen);
    g = (num_each_repl[idx] / sum_num_each_repl) * ngme_repls[idx]->grad();
  }

// if (debug) std::cout << "g in grad() in ngme class = " << g << std::endl;
  return g;
}

void Ngme::burn_in(int iterations) {
  #pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i=0; i < n_repl; i++) {
    ngme_repls[i]->burn_in(iterations);
  }
}

VectorXd Ngme::get_parameter() {
  VectorXd p = ngme_repls[0]->get_parameter();
if (debug) std::cout << "p in get_parameter() in ngme class = " << p << std::endl;
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
}