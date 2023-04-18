// Notice here:
// for R interface, theta_K is directly the parameter_K of the Opeartor object
// for C interface, theta_K is f(parameter_K) such that it is unbounded

#include "../latent.h"

/*
    AR model:
        parameter_K(0) = alpha
        K = C * alpha + G
*/

// W_size = V_size
// get_K_params, grad_K_params, set_K_params, output

AR::AR(const Rcpp::List& model_list, unsigned long seed, Type type)
    : Latent(model_list, seed),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
    type        (type)
{
    // init K and Q
    K = getK(theta_K);
    if (W_size == V_size) {
        lu_solver_K.init(W_size, 0,0,0);
        lu_solver_K.analyze(K);
    }

    SparseMatrix<double> Q = K.transpose() * K;
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);

    update_each_iter();
}

// wrt. parameter_K (bounded parameter)
SparseMatrix<double> AR::getK(const VectorXd& theta_K) const {
    if (type==Type::rw) return C + G;

    assert (theta_K.size() == 1);
    double theta;
    if (type==Type::ar)
        theta = th2a(theta_K(0));

    SparseMatrix<double> K = theta * C + G;
    return K;
}

SparseMatrix<double> AR::get_dK(int index, const VectorXd& alpha) const {
    assert(index==0);
    double th = theta_K(0);
    double da  = 2 * (exp(th) / pow(1+exp(th), 2));
    return da * C;
}
