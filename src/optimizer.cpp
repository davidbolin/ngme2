#include <cmath>
#include "optimizer.h"
//#include <tuple>

VectorXd optimizer::sgd( model& m,
                    const VectorXd& x0, 
                    const double stepsize, 
                    const double eps,
                    const bool precondioner) {

    int count = 0;
    VectorXd x = x0;

    bool terminate = false;

    while (!terminate)
    {
        count += 1;
        VectorXd grad = m.gradient(x);
        if (precondioner) {
            MatrixXd cond = m.hessian(x);

            // update x <- x - stepsize * H^-1 * grad(x)
            x = x - stepsize * cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
//            x = x - cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
        } else {
            x = x - stepsize * grad;
        }

        // to-do: criteria of eps
        if (grad.norm() <= pow(10, -6))
            terminate = true;
    }
    return x;
}
