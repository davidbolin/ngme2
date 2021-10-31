// Store parameters for optimization

#ifndef NGME_MODEL_H
#define NGME_MODEL_H

#include <Eigen/Dense>
#include "Block.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Model
{
friend class optimizer;
private:
// current w and V -> generate new sample
    VectorXd Y, w, V;
    

    string noise;
public:
    // Model(initial value)

    // to-do: how to generate preconditioner
    MatrixXd precond();

    // to-do: generate gradient
    VectorXd grad();


    // -----  Only for testing, f(x, y) = 3x^2 + 2y^2 + x + 3y + 5;
    // gradient of f(x, y)
    VectorXd gradient(const VectorXd& x) {
        VectorXd g (2);
        g(0) = 6 * x(0) + 1;
        g(1) = 4 * x(1) + 3;
        return g;
    }

    // hessian of f(x, y)
    MatrixXd hessian(const VectorXd&) {
        MatrixXd H (2,2);
        H << 6, 0,
             0, 4;
        return H;
    };


};


#endif