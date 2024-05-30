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

AR::AR(const Rcpp::List& operator_list):
    Operator(operator_list),
    G       (Rcpp::as<SparseMatrix<double,0,int>> (operator_list["G"])),
    C       (Rcpp::as<SparseMatrix<double,0,int>> (operator_list["C"]))
{}

// wrt. parameter_K (bounded parameter)
void AR::update_K(const VectorXd& theta_K) {
    assert (theta_K.size() == 1);
    
    // G[1, 1] = sqrt(1-rho**2)
    double rho = th2a(theta_K(0));
    
    // G.coeffRef(0, 0) = sqrt(1 - pow(rho, 2));
    K = rho * C + G;
    K.coeffRef(0, 0) = sqrt(1 - pow(rho, 2));
}

void AR::update_dK(const VectorXd& theta_K) {
    assert(theta_K.size() == 1);

    double th = theta_K(0);
    double rho = th2a(th);
    double da = 2 * (exp(th) / pow(1+exp(th), 2));
    dK[0] = da * C;
    dK[0].coeffRef(0, 0) = da * (-rho) / sqrt(1 - pow(rho, 2));
}
