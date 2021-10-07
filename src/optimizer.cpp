#include <cmath>
#include "optimizer.h"

VectorXd optimizer::gd( model& m, 
                    const VectorXd x0, 
                    const double stepsize, 
                    const double eps,
                    const bool precondion) {
    
    VectorXd x = x0;

    bool terminate = false;

    while (!terminate)
    {
        VectorXd grad = m.gradient(x);
        if (precondion) {
            MatrixXd cond = m.precond(x);
            x = x - stepsize * cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
        } else {
            x = x - stepsize * grad;
        }

        // to-do: criteria of eps
        if (grad.norm() <= pow(10, -6))
            terminate = true;
    }
    return x;
}
