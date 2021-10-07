// Store parameters for optimization

#ifndef NGME_MODEL_H
#define NGME_MODEL_H

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class model
{
friend class optimizer;
private:
public:
    // to-do: place-holder
    VectorXd gradient(VectorXd x) {
        VectorXd g (2);
        g(0) = 2 * x(0) * x(0);
        g(1) = 0.1 * x(1) * x(1);
        return g;
    }

    // to-do: place-holder
    MatrixXd precond(VectorXd) {
        MatrixXd m;
        return m;
    };
};


#endif