#include "rgig.h"
#include <chrono>
#include "sample_rGIG.h"
using Eigen::VectorXd;

// [[Rcpp::export]]
Eigen::VectorXd rGIG_cpp(
  const Eigen::VectorXd& p,
	const Eigen::VectorXd& a,
	const Eigen::VectorXd& b,
  unsigned long     seed
) {

  gig sampler;
  if(seed == 0)
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  sampler.seed(seed);

  VectorXd V(p.size());
  for(int i = 0; i < p.size(); i++)
    V[i] = sampler.sample(p[i], a[i], b[i]);
  return V;
}

double rGIG_cpp(double p,
                       double   a,
                       double   b,
                       unsigned long     seed) {

  gig sampler;
  if(seed == 0)
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  sampler.seed(seed);

  double V = sampler.sample(p, a, b);
  return V;
}
