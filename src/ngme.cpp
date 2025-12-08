// implement the Ngme class and rand effect class
#include "ngme.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#pragma omp declare reduction(vec_plus : Eigen::VectorXd : omp_out += omp_in)  \
    initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))
#pragma omp declare reduction(mat_plus : Eigen::MatrixXd : omp_out += omp_in)  \
    initializer(omp_priv =                                                     \
                    Eigen::MatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))
#endif

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List &R_ngme, unsigned long seed, int sampling_strategy,
           int num_threads_repl, double start_sd)
    : n_params(Rcpp::as<int>(R_ngme["n_params"])),
      sampling_strategy(sampling_strategy), num_threads_repl(num_threads_repl),
      gen(seed), debug(false) {
  Rcpp::List control_ngme = R_ngme["control_ngme"];
  debug = Rcpp::as<bool>(control_ngme["debug"]);
  stepsize =
      VectorXd::Constant(n_params, Rcpp::as<double>(control_ngme["stepsize"]));

  Rcpp::List list_replicates = Rcpp::as<Rcpp::List>(R_ngme["replicates"]);
  n_repl = list_replicates.size();

  // build for each replicates
  for (int i = 0; i < n_repl; i++) {
    Rcpp::List ngme_repl = Rcpp::as<Rcpp::List>(list_replicates[i]);
    ngme_repls.push_back(std::make_shared<BlockModel>(ngme_repl, seed + i));
    num_each_repl.push_back(ngme_repls[i]->get_n_obs());
  }
  sum_num_each_repl =
      std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);

  // Init the random number generator
  weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(),
                                                     num_each_repl.end());

  VectorXd p = ngme_repls[0]->get_parameter();
  std::normal_distribution<double> distribution(0.0, start_sd);
  for (int i = 0; i < p.size(); ++i) {
    p[i] += distribution(gen);
  }
  ngme_repls[0]->set_parameter_and_update(p, true);
}

void Ngme::compute(bool with_precond, double eps) {
  // Aggregate over replicates
  last_grad_ = VectorXd::Zero(n_params);
  // If with_precond=true, also compute and cache preconditioner for current
  // strategy
  if (with_precond) {
    last_precond_ = MatrixXd::Zero(n_params, n_params);
    precond_valid_ = false;
  }
  if (sampling_strategy == Strategy::all) {
    for (int i = 0; i < n_repl; i++) {
      ngme_repls[i]->compute_grad_and_hessian(with_precond, eps);
      last_grad_ += ngme_repls[i]->get_gradient();
      if (with_precond)
        last_precond_ += ngme_repls[i]->get_preconditioner();
    }
  } else { // ws
    int idx = weighted_sampler(gen);
    ngme_repls[idx]->compute_grad_and_hessian(with_precond, eps);
    last_grad_ = ngme_repls[idx]->get_gradient();
    if (with_precond)
      last_precond_ = ngme_repls[idx]->get_preconditioner();
  }
  grad_valid_ = true;
  if (with_precond) {
    precond_valid_ = true;
  }
  mark_computed();
}

MatrixXd Ngme::precond() {
  if (!grad_valid_)
    throw std::runtime_error("precond() called before compute()");
  return -last_precond_;
}

VectorXd Ngme::grad() {
  if (!grad_valid_)
    throw std::runtime_error("grad() called before compute()");
  return -last_grad_;
}

void Ngme::burn_in(int iterations) {
#pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i = 0; i < n_repl; i++) {
    ngme_repls[i]->burn_in(iterations);
  }
}

VectorXd Ngme::get_parameter() {
  VectorXd p = ngme_repls[0]->get_parameter();
  return p;
}

void Ngme::set_parameter_and_update(const VectorXd &p, bool with_precond) {
  // set the same parameter for all latent
  // std::chrono::steady_clock::time_point begin =
  // std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i = 0; i < n_repl; i++) {
    ngme_repls[i]->set_parameter_and_update(p, with_precond);
  }

  // set the different parameter for each random effect
  if (debug)
    std::cout << "set_parameter() in ngme class" << std::endl;
  // Invalidate caches
  grad_valid_ = false;
  precond_valid_ = false;
  reset_computed();
}

bool Ngme::any_latent_nu_at_lower_bound() const {
  for (const auto &repl : ngme_repls) {
    if (repl->any_latent_nu_at_lower_bound()) {
      return true;
    }
  }
  return false;
}

std::vector<std::string> Ngme::latent_nu_lower_bound_summaries() const {
  std::vector<std::string> summaries;
  for (size_t i = 0; i < ngme_repls.size(); ++i) {
    auto repl_summaries = ngme_repls[i]->latent_nu_lower_bound_summaries();
    for (const auto &s : repl_summaries) {
      if (n_repl > 1) {
        std::ostringstream oss;
        oss << "replicate " << i + 1 << ": " << s;
        summaries.push_back(oss.str());
      } else {
        summaries.push_back(s);
      }
    }
  }
  return summaries;
}
