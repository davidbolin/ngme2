// #include <Rcpp.h>

#include "rgig.h"
#include <chrono>
#include <RcppEigen.h>
#include "sample_rGIG.h"

// [[Rcpp::export]]
Eigen::VectorXd rGIG_cpp(Eigen::VectorXd p,
                       Eigen::VectorXd   a,
                       Eigen::VectorXd   b,
                       unsigned long     seed) {

  gig sampler;
  if(seed == 0)
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  sampler.seed(seed);

  Eigen::VectorXd V(p.size());
  for(int i = 0; i < p.size(); i++)
    V[i] = sampler.sample(p[i], a[i], b[i]);
  return V;
}
