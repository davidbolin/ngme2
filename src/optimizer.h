#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include "model.h"
#include <memory>

class Optimizer
{
private:
    bool verbose {false};
    int precond_strategy; // 0 for none, 1 for fast, 2 for full
    double precond_eps;
    int curr_iter;

    int max_iter;
    double max_relative_step;
    double max_absolute_step;
    double eps;
    double stepsize;

    // store the preconditioner
    std::shared_ptr<MatrixXd> precondioner;

    // keep trajs
    std::vector<VectorXd> trajs;
public:
    Optimizer(const Rcpp::List& control_opt);

    // provide model.get_stepsizes()
    // works with procond_grad
     Eigen::VectorXd sgd(
        Model& model,
        double eps,
        int iterations,
        double max_relative_step,
        double max_absolute_step,
        bool compute_precond_each_iter = false
    );

    std::vector<VectorXd> get_trajs() const { return trajs; }

    // getter and setter of preconditioner
    MatrixXd get_precondioner() const { return *precondioner; }
    void set_precondioner(const MatrixXd& precondioner) { this->precondioner = std::make_shared<MatrixXd>(precondioner); }
};

#endif
