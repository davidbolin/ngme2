#ifndef NGME_PRIOR_H
#define NGME_PRIOR_H

#include <Eigen/Dense>
#include <random>
#include <string>
#include <cmath>
#include "sample_rGIG.h"

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

using Eigen::VectorXd;
using std::string;

class PriorUtil {
public:
  static double log_dens(
    const string& type,
    const VectorXd& params,
    double value
  );

  static double d_log_dens(
    const string& type,
    const VectorXd& params,
    double value
  );
};

#endif