#ifndef NGME_OPT_H
#define NGME_OPT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include "model.h"

class Optimizer
{
private:
    bool verbose {false};

    int max_iter;
    double max_relative_step;
    double max_absolute_step;
    double eps;
    double stepsize;
    bool precondioner;
    int curr_iter;

    // keep trajs
    std::vector<VectorXd> trajs;
public:
    Optimizer(const Rcpp::List& control_opt)
        : verbose(control_opt["verbose"]),
        curr_iter(0)
        // max_iter(control_opt["max_iter"]),
        // max_relative_step(control_opt["max_relative_step"]),
        // max_absolute_step(control_opt["max_absolute_step"]),
        // eps(control_opt["eps"]),
        // stepsize(control_opt["stepsize"]),
        // precondioner(control_opt["precondioner"]),
        // iterations(control_opt["iterations"])
    {}

    Rcpp::List sgd(Model& model,
                double stepsize,
                double eps,
                bool precondioner,
                int iterations);

    // provide model.get_stepsizes()
    // works with procond_grad
     Eigen::VectorXd sgd(
            Model& model,
            double eps,
            int iterations,
            double max_relative_step,
            double max_absolute_step);

    std::vector<VectorXd> get_trajs() const { return trajs; }
};

#endif
