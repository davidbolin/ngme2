#ifndef NGME_SAMPLE_RGIG_H
#define NGME_SAMPLE_RGIG_H

#include <RcppEigen.h>
using Eigen::VectorXd;
using Eigen::Ref;

Eigen::VectorXd rGIG_cpp(
	const Eigen::VectorXd&,
	const Eigen::VectorXd&,
	const Eigen::VectorXd&,
	unsigned long=0
);

double rGIG_cpp(double, double, double, unsigned long=0);

#endif