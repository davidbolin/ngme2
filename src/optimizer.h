#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <vector>

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "model.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

using std::vector;

class Optimizer
{
private:
    int iterations;
public:
    // terminate criteria
    // regularizer -> proximal operator wrt regularization
    Rcpp::List sgd(Model& model,
                double stepsize, 
                double eps,
                bool precondioner,
                int iterations);

};

#endif
