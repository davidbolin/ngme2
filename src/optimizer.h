#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "model.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

class Optimizer
{
private:
public:
    // to-do: terminate criteria
    // to-do: regularizer? -> proximal operator wrt regularization
    VectorXd sgd(Model& model,
                const double stepsize, 
                const double eps,
                const bool precondioner);
};

#endif
