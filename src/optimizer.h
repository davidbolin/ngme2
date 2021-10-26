#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "model.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

class optimizer
{
private:
public:
    // to-do: terminate criteria
    // to-do: regularizer? -> proximal operator wrt regularization
    VectorXd sgd(Model& m1,
                const VectorXd& x0, 
                const double stepsize, 
                const double eps,
                const bool precondioner);



};

#endif
