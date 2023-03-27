// implement the Ngme class and rand effect calss
#include "ngme.h"

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List& list_ngmes, unsigned long seed, int sampling_strategy) :
  n_repl(list_ngmes.size()),
  sampling_strategy(sampling_strategy),
  gen(seed)
{
  for (int i=0; i < n_repl; i++) {
    Rcpp::List block_model = Rcpp::as<Rcpp::List> (list_ngmes[i]);
    ngme_repls.push_back(std::make_unique<BlockModel>(block_model, seed));
    num_each_repl.push_back(ngme_repls[i]->get_n_obs());
  }
  sum_num_each_repl = std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);
  n_params = ngme_repls[0]->get_n_params();

  // Init the random number generator
  weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(), num_each_repl.end());
}

VectorXd Ngme::precond_grad() {
    VectorXd g = VectorXd::Zero(n_params);
    // weighted averge over all replicates
    if (sampling_strategy == Strategy::all) {
      for (int i=0; i < n_repl; i++) {
        g +=  (num_each_repl[i] / sum_num_each_repl) * ngme_repls[i]->precond_grad();
// std::cout << "sum_num_each_repl = " << sum_num_each_repl << std::endl;
      }
// std::cout << "using weighted average" << std::endl;
    } else if (sampling_strategy == Strategy::ws) {
      // importance sampling (IS) for each replicate
      int idx = weighted_sampler(gen);
      g =  (num_each_repl[idx] / sum_num_each_repl) * ngme_repls[idx]->precond_grad();
    } else {
      Rcpp::stop("Unknown sampling strategy!");
    }

    return g;
  }