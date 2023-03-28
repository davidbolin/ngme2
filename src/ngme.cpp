// implement the Ngme class and rand effect calss
#include "ngme.h"

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List& R_ngme, unsigned long seed, int sampling_strategy) :
  n_params          (Rcpp::as<int> (R_ngme["n_params"])),
  n_re_params       (Rcpp::as<int> (R_ngme["n_re_params"])),
  sampling_strategy (sampling_strategy),
  gen               (seed),
  debug             (false)
{
  Rcpp::List control_ngme = R_ngme["control_ngme"];
    debug = Rcpp::as<bool> (control_ngme["debug"]);

  Rcpp::List list_replicates = Rcpp::as<Rcpp::List> (R_ngme["replicates"]);
  n_repl = list_replicates.size();

  for (int i=0; i < n_repl; i++) {
    Rcpp::List ngme_repl = Rcpp::as<Rcpp::List> (list_replicates[i]);
    ngme_repls.push_back(std::make_unique<BlockModel>(ngme_repl, seed));
    num_each_repl.push_back(ngme_repls[i]->get_n_obs());
  }
  sum_num_each_repl = std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);

  // Init the random number generator
  weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(), num_each_repl.end());
}

VectorXd Ngme::precond_grad() {
  VectorXd g = VectorXd::Zero(n_params);
  // --- gradient except for random effects ---
  VectorXd g_except_reff = VectorXd::Zero(n_params - n_re_params * n_repl);
  // weighted averge over all replicates
  if (sampling_strategy == Strategy::all) {
    for (int i=0; i < n_repl; i++) {
      g_except_reff +=  (num_each_repl[i] / sum_num_each_repl) * ngme_repls[i]->precond_grad_no_reff();
    }
  } else if (sampling_strategy == Strategy::ws) {
    // weighted sampling (WS) for each replicate
    int idx = weighted_sampler(gen);
    g_except_reff = (num_each_repl[idx] / sum_num_each_repl) * ngme_repls[idx]->precond_grad_no_reff();
  }

  // --- gradient for random effects ---
  VectorXd g_reff = VectorXd::Zero(n_re_params * n_repl);
  for (int i=0; i < n_repl; i++) {
    g_reff.segment(i * n_re_params, n_re_params) = ngme_repls[i]->precond_grad_reff();
  }

if (debug) std::cout << "g in precond_grad() in ngme class = " << g << std::endl;
  // concat g_except_reff and g_reff
  g << g_except_reff, g_reff;
  return g;
}

MatrixXd Ngme::precond() const {
  MatrixXd precond = MatrixXd::Zero(n_params, n_params);
  for (int i=0; i < n_repl; i++) {
    precond += ngme_repls[i]->precond();
  }
  return precond;
}

VectorXd Ngme::grad() {
  // Not implemented yet
  VectorXd g = VectorXd::Zero(n_params);
  // for (int i=0; i < n_repl; i++) {
  //   g += ngme_repls[i]->grad();
  // }
  return g;
}

VectorXd Ngme::get_parameter() const {
  VectorXd p = VectorXd::Zero(n_params);

  VectorXd p_except_reff = ngme_repls[0]->get_parameter_no_reff();
  VectorXd p_reff = VectorXd::Zero(n_re_params * n_repl);
  for (int i=0; i < n_repl; i++) {
    p_reff.segment(i * n_re_params, n_re_params) = ngme_repls[i]->get_parameter_reff();
  }
  p << p_except_reff, p_reff;

if (debug) std::cout << "p in get_parameter() in ngme class = " << p << std::endl;
  return p;
}

void Ngme::set_parameter(const VectorXd& p) {
  // set the same parameter for all latent
  VectorXd p_except_reff = p.segment(0, n_params - n_re_params * n_repl);
  for (int i=0; i < n_repl; i++)
    ngme_repls[i]->set_parameter_no_reff(p_except_reff);

  // set the different parameter for each random effect
  int idx = n_params - n_re_params * n_repl;
  for (int i=0; i < n_repl; i++) {
    ngme_repls[i]->set_parameter_reff(p.segment(idx, n_re_params));
    idx += n_re_params;
  }
  if (debug) std::cout << "set_parameter() in ngme class" << std::endl;
}