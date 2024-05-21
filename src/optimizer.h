#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include "ngme.h"
#include <memory>

class Ngme_optimizer
{
private:
    std::shared_ptr<Ngme> model;

    bool verbose {false};
    int precond_strategy; // 0 for none, 1 for fast, 2 for full
    double numerical_eps;
    int curr_iter;

    int max_iter;
    double max_relative_step;
    double max_absolute_step;
    double stepsize;

    std::string method;
    VectorXd sgd_parameters;

    // variable for different sgd methods
    double beta1 {0.1}, beta2 {0.99}, eps_hat {1e-8}, lambda {0.1};
    VectorXd m, v; // momentum, velocity

    // store the preconditioner
    MatrixXd preconditioner;
    
    // bfgs: H approximates the inverse of Hessian 
    MatrixXd H;
    VectorXd grad, prev_grad, x, prev_x;

    // keep trajs
    std::vector<VectorXd> trajs;
public:
    Ngme_optimizer(
        const Rcpp::List& control_opt,
        std::shared_ptr<Ngme> model
    );

    // provide model.get_stepsizes()
    // works with procond_grad
     Eigen::VectorXd sgd(
        double eps,
        int iterations,
        double max_relative_step,
        double max_absolute_step,
        bool compute_precond_each_iter = false
    );

    std::vector<VectorXd> get_trajs() const { return trajs; }

    // getter and setter of preconditioner
    MatrixXd get_preconditioner() const { return preconditioner; }
    void set_preconditioner(const MatrixXd& preconditioner)
        { this->preconditioner = preconditioner; }
    
    double line_search(
        const VectorXd& x, 
        const VectorXd& p, 
        double phi_0,
        double phi_prime_0,
        double c1 = 1e-4, 
        double c2 = 0.9, 
        double alpha_max = 1.0
    );

    double zoom(
        const VectorXd& x, 
        const VectorXd& p, 
        double alpha_lo, 
        double alpha_hi, 
        double c1, double c2,
        double phi_0,
        double phi_prime_0
    );

    double log_likelihood(const VectorXd& x) { 
        model->set_parameter(x); 
        return model->log_likelihood(); 
    }
    
    // compute numerical gradient using log_likelihood
    VectorXd numerical_grad(const VectorXd& x) {
        VectorXd grad(x.size());
        for (int i=0; i < x.size(); i++) {
            VectorXd x_plus = x;
            x_plus(i) += numerical_eps;
            VectorXd x_minus = x;
            x_minus(i) -= numerical_eps;
            grad(i) = (log_likelihood(x_plus) - log_likelihood(x_minus)) / (2 * numerical_eps);
        }
        return grad;
    }

    // The derivative of log_likelihood in the direction of p
    double directional_derivative(const VectorXd& x, const VectorXd& p) {
        return numerical_grad(x).dot(p);
    }
};

#endif
