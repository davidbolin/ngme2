#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include "model.h"

class Optimizer
{
private:
    int iterations;
public:
    Rcpp::List sgd(Model& model,
                double stepsize,
                double eps,
                bool precondioner,
                int iterations);

    // provide model.get_stepsizes()
     Eigen::VectorXd sgd(
            Model& model,
            double eps,
            int iterations);
};

#endif
