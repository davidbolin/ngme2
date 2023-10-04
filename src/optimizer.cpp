#include <vector>
#include <Eigen/Dense>

#include "include/timer.h"
#include "optimizer.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Ngme_optimizer::Ngme_optimizer(
    const Rcpp::List& control_opt,
    std::shared_ptr<Ngme> ngme
) : model(ngme),
    verbose(control_opt["verbose"]),
    precond_strategy(control_opt["precond_strategy"]),
    precond_eps(control_opt["precond_eps"]),
    curr_iter(0),

    method(Rcpp::as<std::string>(control_opt["sgd_method"])),
    m(VectorXd::Zero(ngme->get_n_params())),
    v(VectorXd::Zero(ngme->get_n_params())),
    preconditioner(ngme->precond(0, precond_eps))
{
    if (method == "adam") {
        sgd_parameters = (Rcpp::as<VectorXd>(control_opt["sgd_parameters"]));
        beta1 = sgd_parameters(0);
        beta2 = sgd_parameters(1);
        eps_hat = sgd_parameters(2);
    } else if (method =="momentum") {
        sgd_parameters = (Rcpp::as<VectorXd>(control_opt["sgd_parameters"]));
        beta1 = sgd_parameters(0);
        beta2 = sgd_parameters(1);
    }
}

// x <- x - model->stepsize() * model->grad()
// return the parameter after sgd
VectorXd Ngme_optimizer::sgd(
    double eps,
    int iterations,
    double max_relative_step, // comparing to x itself
    double max_absolute_step,
    bool compute_precond_each_iter
) {
    int var_reduce_iter = 5000;
    double reduce_power = 0.2; // 0-1, the bigger, the stronger

    // initialize preconditioner
    VectorXd x = model->get_parameter();

// auto timer_grad = std::chrono::steady_clock::now();
    for (int i = 0; i < iterations; i++) {
        trajs.push_back(x);
        VectorXd grad = model->grad();
        if (compute_precond_each_iter) {
            preconditioner = model->precond(precond_strategy, precond_eps);
        }
        grad = preconditioner.llt().solve(grad);

        // which SGD step
        // default: one step = stepsize * H^-1 * grad
        VectorXd one_step;
        if (method == "adam") {
            m = beta1 * m + (1 - beta1) * grad;
            // biased !
            v = beta2 * v + (1 - beta2) * grad.cwiseProduct(grad);
            // VectorXd diag_precond = preconditioner.diagonal();
            // v = beta2 * v + (1 - beta2) * diag_precond;
            VectorXd m_hat = m / (1 - pow(beta1, curr_iter+1));
            VectorXd v_hat = v / (1 - pow(beta2, curr_iter+1));

            one_step = model->get_stepsizes().cwiseProduct(m_hat.cwiseQuotient(v_hat.cwiseSqrt() + VectorXd::Constant(v_hat.size(), eps_hat)));
        } else if (method == "momentum") {
            m = beta1 * m + beta2 * grad;
            one_step = model->get_stepsizes().cwiseProduct(m);
        } else {
            one_step = model->get_stepsizes().cwiseProduct(grad);
        }

/* test NAN
// if (std::isnan(one_step(one_step.size()-1))) {
//     std::cout << "one_step ISNAN = " << one_step << std::endl;
//     return x;
// }
*/

// std::cout << "get gradient (ms): " << since(timer_grad).count() << std::endl;
        // restrict one_step by |one_step(i)| / |x(i)| < rela_step
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
std::cout << "one step = " << one_step << std::endl;
std::cout << "iteration = : " << curr_iter+1 << std::endl;
std::cout << "parameter = : " << x << std::endl;
}

        model->set_parameter(x);
        curr_iter += 1;
    }

    // update preconditioner if not computed
    if (!compute_precond_each_iter)
        preconditioner = model->precond(precond_strategy, precond_eps);

    return x;
}