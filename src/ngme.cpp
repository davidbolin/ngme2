// implement the Ngme class and rand effect calss
#include "ngme.h"

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List& R_ngme, unsigned long seed, int sampling_strategy) :
  n_params          (Rcpp::as<int> (R_ngme["n_params"])),
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

MatrixXd Ngme::precond() const {
  MatrixXd precond = MatrixXd::Zero(n_params, n_params);
  for (int i=0; i < n_repl; i++) {
    precond += ngme_repls[i]->precond();
  }

  return precond / n_repl;
}

VectorXd Ngme::grad() {
  VectorXd g = VectorXd::Zero(n_params);
  // weighted averge over all replicates
  if (sampling_strategy == Strategy::all) {
    for (int i=0; i < n_repl; i++) {
// std::cout << "W repl i " << i << " =" << ngme_repls[i]->getW() << std::endl;
      g +=  (num_each_repl[i] / sum_num_each_repl) * ngme_repls[i]->grad() / n_repl;
    }
  } else if (sampling_strategy == Strategy::ws) {
    // weighted sampling (WS) for each replicate
    int idx = weighted_sampler(gen);
    g = (num_each_repl[idx] / sum_num_each_repl) * ngme_repls[idx]->grad();
  }

if (debug) std::cout << "g in grad() in ngme class = " << g << std::endl;
  return g;
}


VectorXd Ngme::get_parameter() {
  VectorXd p = ngme_repls[0]->get_parameter();
if (debug) std::cout << "p in get_parameter() in ngme class = " << p << std::endl;
  return p;
}

void Ngme::set_parameter(const VectorXd& p) {
  // set the same parameter for all latent
  for (int i=0; i < n_repl; i++)
    ngme_repls[i]->set_parameter(p);

  // set the different parameter for each random effect
  if (debug) std::cout << "set_parameter() in ngme class" << std::endl;
}