#include <Eigen/Dense>
#include <vector>

#include "include/timer.h"
#include "optimizer.h"
#include <sstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Ngme_optimizer::Ngme_optimizer(const Rcpp::List &control_opt,
                               std::shared_ptr<Ngme> ngme)
    : model(ngme), verbose(control_opt["verbose"]),
      numerical_eps(control_opt["numerical_eps"]), curr_iter(0),

      method(Rcpp::as<std::string>(control_opt["sgd_method"])),
      m(VectorXd::Zero(ngme->get_n_params())),
      v(VectorXd::Zero(ngme->get_n_params())),
      preconditioner(
          MatrixXd::Identity(ngme->get_n_params(), ngme->get_n_params())),
      grad(VectorXd::Zero(ngme->get_n_params())), x(ngme->get_parameter()) {
  compute_precond = (method == "precond_sgd");
  if (method != "precond_sgd" && method != "bfgs") {
    // Parse optimizer-specific parameters if any (none for vanilla sgd)
    sgd_parameters = (Rcpp::as<VectorXd>(control_opt["sgd_parameters"]));
  }

  if (method == "momentum") {
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
    line_search_method = Rcpp::as<std::string>(control_opt["line_search"]);
  } else if (method == "sgd") {
    // vanilla sgd requires no extra initialization
  } else {
    // precond_sgd
  }

  // Stepsize decay configuration
  if (control_opt.containsElementNamed("stepsize_decay")) {
    std::string decay =
        Rcpp::as<std::string>(control_opt["stepsize_decay"]);
    stepsize_decay_enabled = (decay == "grad_norm_plateau");
  }
  if (stepsize_decay_enabled) {
    stepsize_decay_min_stepsize =
        Rcpp::as<double>(control_opt["stepsize_decay_min_stepsize"]);
    stepsize_decay_base_stepsize = model->get_stepsizes().mean();
    // Disable for line-search methods
    if (method == "bfgs") {
      stepsize_decay_enabled = false;
    }
  }

  // Do not initialize preconditioner here; build it on-demand inside sgd loop
  x = model->get_parameter();
  prev_x = x;

  model->set_parameter_and_update(x, compute_precond);
}

void Ngme_optimizer::log_verbose_message(const std::string &msg) const {
#ifdef _OPENMP
#pragma omp critical(ngme_verbose_print)
  { Rcpp::Rcout << msg; }
#else
  Rcpp::Rcout << msg;
#endif
}

// x <- x - model->stepsize() * model->grad()
// return the parameter after sgd
VectorXd Ngme_optimizer::sgd(double eps, int iterations,
                             double max_relative_step, // comparing to x itself
                             double max_absolute_step,
                             bool compute_precond_each_iter) {
  int var_reduce_iter = 5000;
  double reduce_power = 0.2; // 0-1, the bigger, the stronger

  // Initialize for adaptive gradient descent
  double avg_stepsize = model->get_stepsizes().mean();
  double theta = 1e9, prev_avg_stepsize = avg_stepsize;

  // auto timer_grad = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; i++) {
    trajs.push_back(x);

    if (method != "bfgs") {
      // Unified compute → then get
      // compute gradients only
      model->compute(compute_precond, numerical_eps);
      grad = model->grad();
    } else {
      grad = numerical_grad(x);
    }

    last_grad_norm = grad.norm();

    // Pflug diagnostic
    if (pflug_conv_check && curr_iter > 0) {
      double inner_prod = grad.dot(prev_grad);
      pflug_sum += inner_prod;
      max_pflug_sum = std::max(max_pflug_sum, pflug_sum);
    }

    // which SGD step
    // default: one step = stepsize * H^-1 * grad
    VectorXd one_step;
    if (method == "adam") {
      m = beta1 * m + (1 - beta1) * grad;
      v = beta2 * v + (1 - beta2) * grad.cwiseProduct(grad);
      VectorXd m_hat = m / (1 - pow(beta1, curr_iter + 1));
      VectorXd v_hat = v / (1 - pow(beta2, curr_iter + 1));
      one_step = model->get_stepsizes().cwiseProduct(m_hat.cwiseQuotient(
          v_hat.cwiseSqrt() + VectorXd::Constant(v_hat.size(), eps_hat)));
    } else if (method == "adamW") {
      m = beta1 * m + (1 - beta1) * grad;
      v = beta2 * v + (1 - beta2) * grad.cwiseProduct(grad);
      VectorXd m_hat = m / (1 - pow(beta1, curr_iter + 1));
      VectorXd v_hat = v / (1 - pow(beta2, curr_iter + 1));
      one_step = model->get_stepsizes().cwiseProduct(
          lambda * x + // extra term for generalization
          m_hat.cwiseQuotient(v_hat.cwiseSqrt() +
                              VectorXd::Constant(v_hat.size(), eps_hat)));
    } else if (method == "momentum") {
      m = beta1 * m + beta2 * grad;
      one_step = model->get_stepsizes().cwiseProduct(m);
    } else if (method == "adagrad") {
      //  v_t = v_{t-1} + g_t^2
      //  x_{t+1} = x_t - stepsize * g_t / (sqrt(v_t) + epsilon)
      v = v + grad.cwiseProduct(grad);
      one_step = model->get_stepsizes().cwiseProduct(grad.cwiseQuotient(
          v.cwiseSqrt() + VectorXd::Constant(v.size(), eps_hat)));
    } else if (method == "rmsprop") {
      //  v_t = beta1 * v_{t-1} + (1-beta1) * g_t^2
      //  x_{t+1} = x_t - stepsize * g_t / (sqrt(v_t) + epsilon)
      v = beta1 * v + (1 - beta1) * grad.cwiseProduct(grad);
      one_step = model->get_stepsizes().cwiseProduct(grad.cwiseQuotient(
          v.cwiseSqrt() + VectorXd::Constant(v.size(), eps_hat)));
    } else if (method == "adaptive_gd") {
      // adaptive gradient descent
      // The update rule for adaptive gradient descent is:
      // \deqn{\lambda_k = \min(\sqrt{1 + \theta_{k-1}} \lambda_{k-1},
      // \frac{||x_k - x_{k-1}||}{2 ||\nabla f(x_k) - \nabla f(x_{k-1})||} )}
      // \deqn{x_{k+1} = x_k - \lambda_k \nabla f(x_k)}
      // \deqn{\theta_k = \lambda_k / \lambda_{k-1}}
      if (i == 0) {
        one_step = avg_stepsize * grad;
      } else {
        avg_stepsize =
            std::min(std::sqrt(1 + theta) * prev_avg_stepsize,
                     0.5 * (x - prev_x).norm() / (grad - prev_grad).norm());
        theta = avg_stepsize / prev_avg_stepsize;
        prev_grad = grad;
        one_step = avg_stepsize * grad;
      }
      // std::cout << "avg_stepsize = " << avg_stepsize << std::endl;
      // std::cout << "theta = " << theta << std::endl;
    } else if (method == "bfgs" && curr_iter > 0) {
      // BFGS - deterministic for Gaussian model
      VectorXd s = x - prev_x;
      VectorXd y = grad - prev_grad;
      double rho = 1.0 / y.dot(s);
      H = (MatrixXd::Identity(H.rows(), H.cols()) - rho * s * y.transpose()) *
              H *
              (MatrixXd::Identity(H.rows(), H.cols()) -
               rho * y * s.transpose()) +
          rho * s * s.transpose();
      // compute direction
      one_step = -H * grad;

      // make sure the direction is descent
      if (grad.dot(one_step) > 0) {
        one_step = -one_step;
      }

      double curr_loglik = model->log_likelihood();
      double alpha = line_search(line_search_method, x, one_step, curr_loglik,
                                 grad.dot(one_step),
                                 1e-4, // c1
                                 0.9,  // c2
                                 1.1   // alpha_max
      );

      // note here feed the opposite direction
      one_step = -alpha * one_step;
      prev_x = x;
      prev_grad = grad;
    } else if (method == "precond_sgd") {
      // Always use full preconditioner for precond_sgd
      preconditioner = model->precond();
      grad = preconditioner.llt().solve(grad);
      one_step = model->get_stepsizes().cwiseProduct(grad);
    } else if (method == "sgd") {
      // Vanilla SGD (no preconditioner)
      one_step = model->get_stepsizes().cwiseProduct(grad);
    }

    if (stepsize_decay_enabled && method != "bfgs") {
      one_step *= stepsize_decay_scale;
    }

    // Test if one_step is NAN
    if (std::isnan(one_step(one_step.size() - 1))) {
      std::ostringstream oss;
      oss << "grad.norm() = " << grad.norm() << '\n';
      oss << " H = " << H << '\n';
      oss << "one_step ISNAN = " << one_step << '\n';
      log_verbose_message(oss.str());
      return x;
    }

    // std::cout << "get gradient (ms): " << since(timer_grad).count() <<
    // std::endl; restrict one_step by |one_step(i)| / |x(i)| < rela_step
    VectorXd rela_max_step = max_relative_step * x.cwiseAbs();
    for (int j = 0; j < one_step.size(); j++) {
      double sign = one_step(j) > 0 ? 1.0 : -1.0;

      // // take limit on relative step
      // if (abs(x(j)) > 1 && abs(one_step(j)) > rela_max_step(j)) {
      //     one_step(j) = sign * rela_max_step(j);
      // }

      // take limit on absolute step
      if (abs(one_step(j)) > max_absolute_step) {
        one_step(j) = sign * max_absolute_step;
      }
    }

    // variance reduction (not enabled)
    int r = curr_iter - var_reduce_iter;
    double tmp = r > 0 ? (1.0 / r) : 1.0;
    // x = x - pow(tmp, reduce_power) * one_step;
    x = x - one_step;

    if (verbose) {
      std::ostringstream oss;
      oss << "iteration = : " << curr_iter + 1 << '\n';
      oss << "grad.norm() = " << grad.norm() << '\n';
      if (stepsize_decay_enabled) {
        oss << "stepsize_scale = " << stepsize_decay_scale << '\n';
        if (stepsize_decay_base_stepsize > 0.0) {
          oss << "effective_stepsize = "
              << (stepsize_decay_scale * stepsize_decay_base_stepsize) << '\n';
        }
      }
      // oss << "parameter = : " << x << '\n';
      // oss << "marginal likelihood := " <<  -model->log_likelihood() << '\n';
      if (pflug_conv_check) {
        oss << "pflug_sum = " << pflug_sum
            << ", max_pflug_sum = " << max_pflug_sum << '\n';
      }
      oss << "---------------------------\n";
      log_verbose_message(oss.str());
    }

    model->set_parameter_and_update(x, compute_precond);
    prev_grad = grad;
    curr_iter += 1;
  }

  return x;
}

// line_search algo for BFGS
// Algorithm 3.5 in Nocedal and Wright
double Ngme_optimizer::line_search_wolfe(const VectorXd &x, const VectorXd &p,
                                         double phi_0, double phi_prime_0,
                                         double c1, double c2,
                                         double alpha_max) {
  // alpha_max = 2;

  double alpha_0 = 0.0;
  // double alpha_i = 0.5 * alpha_max; // Initial step length
  double alpha_i = 1;

  // Initial values
  double phi_alpha_i = phi_0;
  double phi_prime_alpha_i = phi_prime_0;
  double prev_phi_alpha_i = phi_0;
  double prev_alpha_i = 0.0;

  int i = 1;

  while (true) {
    VectorXd x_new = x + alpha_i * p;
    prev_phi_alpha_i = phi_alpha_i;
    phi_alpha_i = log_likelihood(x_new);

    if (phi_alpha_i > phi_0 + c1 * alpha_i * phi_prime_0 ||
        (phi_alpha_i >= prev_phi_alpha_i && i > 1)) {
      return zoom(x, p, alpha_i - 0.5 * alpha_max, alpha_i, c1,
                  c2, // prev_alpha_i, alpha_i
                  phi_0, phi_prime_0);
    }

    phi_prime_alpha_i = directional_derivative(x_new, p);
    if (std::abs(phi_prime_alpha_i) <= -c2 * phi_prime_0) {
      return alpha_i;
    }
    if (phi_prime_alpha_i >= 0) {
      return zoom(x, p, alpha_i, alpha_i - 0.5 * alpha_max, c1, c2, phi_0,
                  phi_prime_0);
    }

    prev_alpha_i = alpha_i;
    alpha_i = 0.5 * (alpha_i + alpha_max); // Update alpha_i
    i++;
    // std::cout << "line search iteration = " << i << std::endl;

    // Termination condition to prevent infinite loop (can be tuned)
    if (i > 50) {
      // std::cout << "warning: line search iteration > 50" << std::endl;
      // std::cout << "alpha_i = " << alpha_i << std::endl;
      return alpha_i;
    }
  }
}

// Zoom function (Algorithm 3.6)
double Ngme_optimizer::zoom(const VectorXd &x, const VectorXd &p,
                            double alpha_lo, double alpha_hi, double c1,
                            double c2, double phi_0, double phi_prime_0) {
  // if (alpha_lo > alpha_hi) {
  //     std::cout << "warning: alpha_lo > alpha_hi" << std::endl;
  // }
  while (true) {
    // Interpolate (using bisection as an example)
    double alpha_j = 0.5 * (alpha_lo + alpha_hi);

    VectorXd x_new = x + alpha_j * p;

    double phi_alpha_j = log_likelihood(x_new);

    if (phi_alpha_j > phi_0 + c1 * alpha_j * phi_prime_0 ||
        phi_alpha_j >= log_likelihood(x + alpha_lo * p)) {
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
      // std::cout << "zoom out" << std::endl;
      return alpha_j;
    }
  }
}

double Ngme_optimizer::line_search_backtracking(const VectorXd &x,
                                                const VectorXd &p, double phi_0,
                                                double phi_prime_0, double c1,
                                                double c2, double alpha_max) {
  double alpha = alpha_max;
  double phi_alpha = phi_0;
  double phi_prime_alpha = phi_prime_0;

  while (phi_alpha > phi_0 + c1 * alpha * phi_prime_0) {
    alpha *= 0.9;
    VectorXd x_new = x + alpha * p;
    phi_alpha = log_likelihood(x_new);
    phi_prime_alpha = directional_derivative(x_new, p);

    if (alpha < 1e-7) {
      // std::cout << "warning: alpha < 1e-7" << std::endl;
      // cannot return 0, H will be NAN
      return alpha;
    }
  }
  return alpha;
}
