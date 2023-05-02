// Notice here:
// for R interface, theta_K is directly the parameter_K of the Opeartor object
// for C interface, theta_K is f(parameter_K) such that it is unbounded

#include "../operator.h"

/*
    AR model:
        parameter_K(0) = alpha
        K = C * alpha + G
*/

// W_size = V_size
// get_K_params, grad_K_params, set_K_params, output

AR::AR(const Rcpp::List& operator_list, Type type):
    Operator(operator_list),
    G       (Rcpp::as<SparseMatrix<double,0,int>> (operator_list["G"])),
    C       (Rcpp::as<SparseMatrix<double,0,int>> (operator_list["C"])),
    type    (type)
{}

// wrt. parameter_K (bounded parameter)
SparseMatrix<double> AR::getK(const VectorXd& theta_K) const {
    if (type==Type::rw) return C + G;

    assert (theta_K.size() == 1);

    if (type==Type::ar) {
        return th2a(theta_K(0)) * C + G;
    } else {
        throw std::invalid_argument("wrong type");
    }
}

SparseMatrix<double> AR::get_dK(int index, const VectorXd& theta_K) const {
    assert(index==0);
    double th = theta_K(0);
    double da = 2 * (exp(th) / pow(1+exp(th), 2));
    return da * C;
}
