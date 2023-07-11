#ifndef NGME_DIFF_H
#define NGME_DIFF_H

#include <Eigen/Dense>
using Eigen::VectorXd, Eigen::MatrixXd;

VectorXd num_g(const VectorXd& v, double (*f)(const VectorXd&), double eps);
MatrixXd num_h(const VectorXd& v, double (*f)(const VectorXd&), double eps);

#endif