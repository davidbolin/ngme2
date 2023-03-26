// implement the Ngme class and rand effect calss
#include "ngme.h"

VectorXd Ngme::precond_grad() {
    VectorXd g = VectorXd::Zero(n_params);
    // weighted averge over all replicates
    if (sampling_strategy == Strategy::all) {
      for (int i=0; i < n_repl; i++) {
        g +=  (num_each_repl[i] / sum_num_each_repl) * blocks[i]->precond_grad();
// std::cout << "sum_num_each_repl = " << sum_num_each_repl << std::endl;
      }
// std::cout << "using weighted average" << std::endl;
    } else if (sampling_strategy == Strategy::ws) {
      // importance sampling (IS) for each replicate
      int idx = weighted_sampler(gen);
      g =  (num_each_repl[idx] / sum_num_each_repl) * blocks[idx]->precond_grad();
    } else {
      Rcpp::stop("Unknown sampling strategy!");
    }

    return g;
  }