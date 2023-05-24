#ifndef NGME_SAMPLE_RGIG_H
#define NGME_SAMPLE_RGIG_H

#include <RcppEigen.h>

Eigen::VectorXd rGIG_cpp(Eigen::VectorXd,
                       	 Eigen::VectorXd,
                       	 Eigen::VectorXd,
                       	 unsigned long=0);

double rGIG_cpp(double, double, double, unsigned long=0);

#endif