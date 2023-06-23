
#include "num_diff.h"

VectorXd num_g(const VectorXd& v, double (*f)(const VectorXd&), double eps) {
	int n = v.size();
	VectorXd grad(n);
	double original_val = f(v);
	for (int i=0; i < n; i++) {
		VectorXd tmp_v = v; tmp_v(i) += eps;
		double tmp_val = f(tmp_v);
		grad(i) = (tmp_val - original_val) / eps;
	}
	return grad;
}

MatrixXd num_h(const VectorXd& v, double (*f)(const VectorXd&), double eps) {
	int n = v.size();
	MatrixXd hessian(n, n);
	double original_val = f(v);
	VectorXd f_v (n);
	// compute f_v = f(v + eps * e_i)
	for (int i=0; i < n; i++) {
		VectorXd tmp_v = v; tmp_v(i) += eps;
		f_v(i) = f(tmp_v);
	}
	// compute H_ij = d2 f / dxi dxj
	for (int i=0; i < n; i++) {
		for (int j=0; j <= i; j++) {
			VectorXd tmp_vij = v; tmp_vij(i) += eps; tmp_vij(j) += eps;
			double f_vij = f(tmp_vij);
			hessian(i, j) = (f_vij - f_v(i) - f_v(j) + original_val) / (eps * eps);
		}
	}
	// fill in the lower triangular part
	for (int i=0; i < n; i++) {
		for (int j=0; j < i; j++) {
			hessian(j, i) = hessian(i, j);
		}
	}
	return hessian;
}