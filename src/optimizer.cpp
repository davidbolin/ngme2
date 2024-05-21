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
    numerical_eps(control_opt["numerical_eps"]),
    curr_iter(0),

    method(Rcpp::as<std::string>(control_opt["sgd_method"])),
    m(VectorXd::Zero(ngme->get_n_params())),
    v(VectorXd::Zero(ngme->get_n_params())),
    preconditioner(ngme->precond(0, numerical_eps)),
    grad(VectorXd::Zero(ngme->get_n_params())),
    x(ngme->get_parameter())
{
    if (method != "precond_sgd" && method != "bfgs") {
        sgd_parameters = (Rcpp::as<VectorXd>(control_opt["sgd_parameters"]));
    }

    if (method =="momentum") {
        beta1 = sgd_parameters(0);
        beta2 = sgd_parameters(1);
    } else if (method == "adagrad") {
        eps_hat = sgd_parameters(0);
    } else if (method == "rmsprop") {
        beta1 = sgd_parameters(0);
        eps_hat = sgd_parameters(1);
    } else if (method == "adam") {
        beta1 = sgd_parameters(0);
        beta2 = sgd_parameters(1);
        eps_hat = sgd_parameters(2);
    } else if (method == "adamW") {
        beta1 = sgd_parameters(0);
        beta2 = sgd_parameters(1);
        lambda = sgd_parameters(2);
        eps_hat = sgd_parameters(3);
    } else if (method == "bfgs") {
        // init for BFGS
        H = MatrixXd::Identity(ngme->get_n_params(), ngme->get_n_params()) * 1e-8;
    } else {
        // precond_sgd
    }

    grad = model->grad();
    prev_grad = grad;
    x = model->get_parameter();
    prev_x = x;
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

// auto timer_grad = std::chrono::steady_clock::now();
    for (int i = 0; i < iterations; i++) {
        trajs.push_back(x);

        if (method != "bfgs") {
            // stochastic gradient descent
            grad = model->grad();
        } else {
            grad = numerical_grad(x);
        }

        // which SGD step
        // default: one step = stepsize * H^-1 * grad
        VectorXd one_step;
        if (method == "adam") {
            m = beta1 * m + (1 - beta1) * grad;
            v = beta2 * v + (1 - beta2) * grad.cwiseProduct(grad);
            VectorXd m_hat = m / (1 - pow(beta1, curr_iter+1));
            VectorXd v_hat = v / (1 - pow(beta2, curr_iter+1));
            one_step = model->get_stepsizes().cwiseProduct(
                m_hat.cwiseQuotient(v_hat.cwiseSqrt() + VectorXd::Constant(v_hat.size(), eps_hat)));
        } else if (method == "adamW") {
            m = beta1 * m + (1 - beta1) * grad;
            v = beta2 * v + (1 - beta2) * grad.cwiseProduct(grad);
            VectorXd m_hat = m / (1 - pow(beta1, curr_iter+1));
            VectorXd v_hat = v / (1 - pow(beta2, curr_iter+1));
            one_step = model->get_stepsizes().cwiseProduct(
                lambda * x + // extra term for generalization
                m_hat.cwiseQuotient(v_hat.cwiseSqrt() + VectorXd::Constant(v_hat.size(), eps_hat)));
        } else if (method == "momentum") {
            m = beta1 * m + beta2 * grad;
            one_step = model->get_stepsizes().cwiseProduct(m);
        } else if (method == "adagrad") {
//  v_t = v_{t-1} + g_t^2
//  x_{t+1} = x_t - stepsize * g_t / (sqrt(v_t) + epsilon)
            v = v + grad.cwiseProduct(grad);
            one_step = model->get_stepsizes().cwiseProduct(grad.cwiseQuotient(v.cwiseSqrt() + VectorXd::Constant(v.size(), eps_hat)));
        } else if (method == "rmsprop") {
//  v_t = beta1 * v_{t-1} + (1-beta1) * g_t^2
//  x_{t+1} = x_t - stepsize * g_t / (sqrt(v_t) + epsilon) 
            v = beta1 * v + (1-beta1) * grad.cwiseProduct(grad);
            one_step = model->get_stepsizes().cwiseProduct(grad.cwiseQuotient(v.cwiseSqrt() + VectorXd::Constant(v.size(), eps_hat)));
        } else if (method == "bfgs" && curr_iter > 0) {
            // BFGS - deterministic for Gaussian model
            VectorXd s = x - prev_x;
            VectorXd y = grad - prev_grad;
            double rho = 1.0 / y.dot(s);
            H = (MatrixXd::Identity(H.rows(), H.cols()) - rho * s * y.transpose()) * H * (MatrixXd::Identity(H.rows(), H.cols()) - rho * y * s.transpose()) + rho * s * s.transpose();
            // compute direction
            one_step = - model->get_stepsizes().cwiseProduct(H * grad);

            double curr_loglik = model->log_likelihood();
            double alpha = line_search(
                x, one_step, curr_loglik, grad.dot(one_step),
                1e-4, // c1
                0.9, // c2
                1.1 // alpha_max
            );
            
            // note here feed the opposite direction
            one_step = - alpha * one_step;
            prev_x = x; prev_grad = grad;
        } else {
            // precond_sgd
            if (compute_precond_each_iter) {
                preconditioner = model->precond(precond_strategy, numerical_eps);
            }
            
            if (precond_strategy > 0) {
                grad = preconditioner.llt().solve(grad);
            }
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
std::cout << "iteration = : " << curr_iter+1 << std::endl;
// std::cout << "one step = " << one_step << std::endl;
// std::cout << "parameter = : " << x << std::endl;
std::cout << "negative marginal likelihood := " <<  model->log_likelihood() << std::endl;
std::cout << "---------------------------" << std::endl; 
}

        model->set_parameter(x);
        curr_iter += 1;
    }

    // update preconditioner if not computed
    if (!compute_precond_each_iter)
        preconditioner = model->precond(precond_strategy, numerical_eps);

    return x;
}

// line_search algo for BFGS
// Algorithm 3.5 in Nocedal and Wright
double Ngme_optimizer::line_search(
    const VectorXd& x, 
    const VectorXd& p, 
    double phi_0,
    double phi_prime_0,
    double c1,   
    double c2,   
    double alpha_max
) {
    double alpha_0 = 0.0;
    // double alpha_i = 0.5 * alpha_max; // Initial step length
    double alpha_i = 1;
    int i = 1;

    while (true) {
        VectorXd x_new = x + alpha_i * p;
        double phi_alpha_i = log_likelihood(x_new);

        if (phi_alpha_i > phi_0 + c1 * alpha_i * phi_prime_0 
        || (phi_alpha_i >= log_likelihood(x + (alpha_i - 0.5 * alpha_max) * p) && i > 1)) {
            return zoom(
                x, p, alpha_i - 0.5 * alpha_max, alpha_i, c1, c2,
                phi_0, phi_prime_0
            );
        }

        double phi_prime_alpha_i = directional_derivative(x_new, p);
        if (std::abs(phi_prime_alpha_i) <= -c2 * phi_prime_0) {
            return alpha_i;
        }
        if (phi_prime_alpha_i >= 0) {
            return zoom(
                x, p, alpha_i, alpha_i - 0.5 * alpha_max, c1, c2,
                phi_0, phi_prime_0
            );
        }

        alpha_i = 0.5 * (alpha_i + alpha_max); // Update alpha_i
        i++;
    }
}


// Zoom function (Algorithm 3.6)
double Ngme_optimizer::zoom(
    const VectorXd& x, 
    const VectorXd& p, 
    double alpha_lo, 
    double alpha_hi, 
    double c1, double c2,
    double phi_0,
    double phi_prime_0
) {
// std::cout << "zoom" << std::endl;
    while (true) {
        // Interpolate (using bisection as an example)
        double alpha_j = 0.5 * (alpha_lo + alpha_hi);

        VectorXd x_new = x + alpha_j * p;
        
        double phi_alpha_j = log_likelihood(x_new);
        
        if (phi_alpha_j > phi_0 + c1 * alpha_j * phi_prime_0 || phi_alpha_j >= log_likelihood(x + alpha_lo * p)) {
            alpha_hi = alpha_j;
        } else {
            double phi_prime_alpha_j = directional_derivative(x_new, p);
            if (std::abs(phi_prime_alpha_j) <= -c2 * phi_prime_0) {
                return alpha_j;
            }
            if (phi_prime_alpha_j * (alpha_hi - alpha_lo) >= 0) {
                alpha_hi = alpha_lo;
            }
            alpha_lo = alpha_j;
        }

        // Termination condition to prevent infinite loop (can be tuned)
        if (std::abs(alpha_hi - alpha_lo) < 1e-8) {
            return alpha_j;
        }
    }
// std::cout << "zoom out" << std::endl;
}