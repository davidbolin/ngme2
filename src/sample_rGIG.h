#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "rgig.h"

Eigen::VectorXd rGIG_cpp(Eigen::VectorXd,
                       	 Eigen::VectorXd,
                       	 Eigen::VectorXd,
                       	 unsigned long=0);