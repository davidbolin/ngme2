// implement the Ngme class and rand effect calss
#include "ngme.h"

VectorXd Ngme::precond_grad() {
    VectorXd g = VectorXd::Zero(n_params);
    // weighted averge over all replicates
    if (sampling_strategy == 1) {
      for (int i=0; i < n_blocks; i++) {
        g += n_obs[i] * blocks[i]->grad() / sum_n_obs;
      }
// std::cout << "using weighted average" << std::endl;
    } else if (sampling_strategy == 2) {
      // importance sampling (IS) for each replicate
      int idx = weighted_sampler(gen);
      g = blocks[idx]->precond_grad() / n_obs[idx];
// std::cout << "using is, idx = " << idx << std::endl;
    } else {
      Rcpp::stop("Unknown sampling strategy!");
    }

    return g;
  }