#include <vector>
#include <Eigen/Dense>

#include "include/timer.h"
#include "optimizer.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Rcpp::List Optimizer::sgd(
    Model& model,
    double stepsize,
    double eps,
    bool precondioner,
    int iterations
) {
    vector<VectorXd> x_traj;
    vector<VectorXd> grad_traj;

    int count = 0;
    VectorXd x = model.get_parameter();

    bool terminate = false;

    while (!terminate)
    {
        count += 1;
// auto timer_grad = std::chrono::steady_clock::now();
        VectorXd grad = model.grad();
// std::cout << "get gradient (ms): " << since(timer_grad).count() << std::endl;

        if (precondioner) {
            MatrixXd cond = model.precond();

            // update x <- x - stepsize * H^-1 * grad(x)
            x = x - stepsize * cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
//            x = x - cond.selfadjointView<Eigen::Upper>().llt().solve(grad);
        } else {
            x = x - stepsize * grad;
        }

        // record x and grad
        x_traj.push_back(x);
        grad_traj.push_back(grad);

        model.set_parameter(x);

        // to-do: criteria of eps
        if ((grad.norm() <= pow(10, -6)) || (count > iterations))
            terminate = true;

    }
    return Rcpp::List::create(Rcpp::Named("grad_traj") = grad_traj,
                              Rcpp::Named("x_traj") = x_traj);
}


// x <- x - model.stepsize() * model.grad()
// return the parameter after sgd
VectorXd Optimizer::sgd(
    Model& model,
    double eps,
    int iterations,
    double max_relative_step, // comparing to x itself
    double max_absolute_step
) {
    // update later
    int var_reduce_iter = 5000;
    double reduce_power = 0.2; // 0-1, the bigger, the stronger

    VectorXd x = model.get_parameter();

    VectorXd grad;

    for (int i = 0; i < iterations; i++) {
        trajs.push_back(x);
// auto timer_grad = std::chrono::steady_clock::now();
        grad = model.precond_grad();
// std::cout << "get gradient (ms): " << since(timer_grad).count() << std::endl;

        // VectorXd stepsizes = model.get_stepsizes();
        // x = x - grad.cwiseProduct(stepsizes);
        // restrict one_step by |one_step(i)| / |x(i)| < rela_step
        VectorXd one_step = grad.cwiseProduct(model.get_stepsizes());

        VectorXd rela_max_step =  max_relative_step * x.cwiseAbs();
        for (int j = 0; j < one_step.size(); j++) {
            double sign = one_step(j) > 0 ? 1.0 : -1.0;

            // take limit on relative step
            if (abs(x(j)) > 1 && abs(one_step(j)) > rela_max_step(j)) {
                one_step(j) = sign * rela_max_step(j);
            }

            // take limit on absolute step
            if (abs(one_step(j)) > max_absolute_step) {
                one_step(j) = sign * max_absolute_step;
            }
        }

        // variance reduction
        int r = curr_iter - var_reduce_iter;
        double tmp = r > 0 ? (1.0/r) : 1.0;

        x = x - pow(tmp, reduce_power) * one_step;

if (verbose) {
std::cout << "iteration = : " << curr_iter << std::endl;
std::cout << "parameter = : " << x << std::endl;
}

        model.set_parameter(x);
        curr_iter += 1;
    }

    return x;
}


        // VectorXd one_step = grad.cwiseProduct(model.get_stepsizes());

        // // restrict one_step by |one_step(i)| / |x(i)| < rela_step
        // // VectorXd ratio = one_step.cwiseAbs().cwiseQuotient(x.cwiseAbs());
        // // for (int j = 0; j < ratio.size() && ratio(j) > max_relative_step; j++) {
        // //     double sign = one_step(j) / abs(one_step(j));
        // //     one_step(j) = sign * max_relative_step * x(j);
        // // }

        // x = x - one_step;
        // model.set_parameter(x);