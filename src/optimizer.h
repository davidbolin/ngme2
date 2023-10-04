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
    double precond_eps;
    int curr_iter;

    int max_iter;
    double max_relative_step;
    double max_absolute_step;
    double eps;
    double stepsize;

    std::string method;
    VectorXd sgd_parameters;

    // ADAM
    double beta1 {0.1}, beta2 {0.99}, eps_hat {1e-8};
    VectorXd m, v; // momentum, velocity

    // store the preconditioner
    MatrixXd preconditioner;

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
};

#endif
