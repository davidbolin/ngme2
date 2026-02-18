#ifndef NGME_OPT_H
#define NGME_OPT_H

#include "ngme.h"
#include <memory>
#include <limits>
#include <random>

class Ngme_optimizer {
private:
  std::shared_ptr<Ngme> model;

  bool verbose{false};
  double numerical_eps;
  int curr_iter;

  int max_iter;
  double max_relative_step;
  double max_absolute_step;
  double stepsize;

  std::string method, line_search_method;
  VectorXd sgd_parameters;

  // variable for different sgd methods
  double beta1{0.1}, beta2{0.99}, eps_hat{1e-8}, lambda{0.1};
  double sgld_temperature{1.0};
  VectorXd m, v; // momentum, velocity

  // store the preconditioner
  MatrixXd preconditioner;
  bool compute_precond{false};

  // bfgs: H approximates the inverse of Hessian
  MatrixXd H;
  VectorXd grad, prev_grad, x, prev_x;

  // keep trajs
  std::vector<VectorXd> trajs;
  bool record_traj{true};

  // Stepsize decay (grad.norm plateau) - scale applied to step
  bool stepsize_decay_enabled{false};
  double stepsize_decay_min_stepsize{0.0};
  double stepsize_decay_scale{1.0};
  double stepsize_decay_base_stepsize{0.0};
  // Stepsize schedule (iteration-based) - e.g. eta_t = eta_0 (t + t0)^(-alpha)
  bool stepsize_schedule_enabled{false};
  std::string stepsize_schedule_method{"constant"};
  double stepsize_schedule_alpha{0.501};
  double stepsize_schedule_t0{1.0};
  int stepsize_schedule_burnin_iter{0};
  double last_grad_norm{0.0};

  // Pflug diagnostic
  bool pflug_conv_check{false};
  double pflug_sum{0.0};
  double max_pflug_sum{0.0};
  std::mt19937 sgld_rng;
  std::normal_distribution<double> sgld_std_normal{0.0, 1.0};

  void log_verbose_message(const std::string &msg) const;

public:
  Ngme_optimizer(const Rcpp::List &control_opt, std::shared_ptr<Ngme> model,
                 unsigned long seed = 0);

  // provide model.get_stepsizes()
  // works with procond_grad
  Eigen::VectorXd sgd(double eps, int iterations, double max_relative_step,
                      double max_absolute_step,
                      bool compute_precond_each_iter = false);

  void set_verbose(bool value) { verbose = value; }
  bool is_verbose() const { return verbose; }

  void set_pflug_conv_check(bool value) { pflug_conv_check = value; }
  bool get_pflug_conv_check() const { return pflug_conv_check; }
  double get_pflug_sum() const { return pflug_sum; }
  double get_max_pflug_sum() const { return max_pflug_sum; }
  double get_last_grad_norm() const { return last_grad_norm; }
  void set_stepsize_decay_scale(double scale) {
    stepsize_decay_scale = scale;
    if (stepsize_decay_min_stepsize > 0.0 &&
        stepsize_decay_base_stepsize > 0.0) {
      double min_scale =
          stepsize_decay_min_stepsize / stepsize_decay_base_stepsize;
      if (min_scale > 1.0) {
        min_scale = 1.0;
      }
      if (stepsize_decay_scale < min_scale) {
        stepsize_decay_scale = min_scale;
      }
    }
  }

  void set_store_traj(bool value) {
    record_traj = value;
    if (!record_traj)
      trajs.clear();
  }
  bool stores_traj() const { return record_traj; }
  std::vector<VectorXd> get_trajs() const { return trajs; }

  void record_current_state() {
    if (record_traj)
      trajs.push_back(x);
  }

  // getter and setter of preconditioner
  MatrixXd get_preconditioner() const { return preconditioner; }
  void set_preconditioner(const MatrixXd &preconditioner) {
    this->preconditioner = preconditioner;
  }

  double line_search_wolfe(const VectorXd &x, const VectorXd &p, double phi_0,
                           double phi_prime_0, double c1 = 1e-4,
                           double c2 = 0.9, double alpha_max = 1.0);

  double line_search_backtracking(const VectorXd &x, const VectorXd &p,
                                  double phi_0, double phi_prime_0,
                                  double c1 = 1e-4, double c2 = 0.9,
                                  double alpha_max = 1.0);

  double zoom(const VectorXd &x, const VectorXd &p, double alpha_lo,
              double alpha_hi, double c1, double c2, double phi_0,
              double phi_prime_0);

  double log_likelihood(const VectorXd &x) {
    model->set_parameter_and_update(x, false);
    return model->log_likelihood();
  }

  // compute numerical gradient using log_likelihood
  VectorXd numerical_grad(const VectorXd &x) {
    VectorXd grad(x.size());
    for (int i = 0; i < x.size(); i++) {
      VectorXd x_plus = x;
      x_plus(i) += numerical_eps;
      VectorXd x_minus = x;
      x_minus(i) -= numerical_eps;
      grad(i) = (log_likelihood(x_plus) - log_likelihood(x_minus)) /
                (2 * numerical_eps);
    }
    return grad;
  }

  // The derivative of log_likelihood in the direction of p
  double directional_derivative(const VectorXd &x, const VectorXd &p) {
    return numerical_grad(x).dot(p);
  }

  double line_search(string method, const VectorXd &x, const VectorXd &p,
                     double phi_0, double phi_prime_0, double c1 = 1e-4,
                     double c2 = 0.9, double alpha_max = 1.0) {
    if (method == "wolfe") {
      return line_search_wolfe(x, p, phi_0, phi_prime_0, c1, c2, alpha_max);
    } else if (method == "backtracking") {
      return line_search_backtracking(x, p, phi_0, phi_prime_0, c1, c2,
                                      alpha_max);
    } else {
      Rcpp::stop("line search method not implemented");
    }
  }
};

#endif
