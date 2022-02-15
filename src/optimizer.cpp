#include "optimizer.h"

VectorXd 
Optimizer::sgd(Model& model,
                        const double stepsize, 
                        const double eps,
                        const bool precondioner,
                        const int iterations) {

    int count = 0;
    VectorXd x = model.get_parameter();

    bool terminate = false;

    while (!terminate)
    {
        count += 1;
        
        VectorXd grad = model.grad();
        
        if (precondioner) {
            MatrixXd cond = model.precond();

            // update x <- x - stepsize * H^-1 * grad(x)
            x = x - stepsize * cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
//            x = x - cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
        } else {
            x = x - stepsize * grad;
        }

        model.set_parameter(x);

        // to-do: criteria of eps
        if ((grad.norm() <= pow(10, -6)) || (count > iterations))
            terminate = true;

    }
    return x;
}