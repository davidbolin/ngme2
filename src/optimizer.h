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
    VectorXd gd(model& m1, 
                const VectorXd x0, 
                const double stepsize, 
                const double eps,
                const bool precondion);

};


#endif
