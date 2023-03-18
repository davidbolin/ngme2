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
if (debug) std::cout << "Begin Constructor of AR1" << std::endl;

    // Init K and Q
    K = getK(theta_K);
    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);

    SparseMatrix<double> Q = K.transpose() * K;

    // watch out!
    if (W_size == V_size) {
        lu_solver_K.init(W_size, 0,0,0);
        lu_solver_K.analyze(K);
        // trace == 0
        // compute_trace();
    }

    // Init Q
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);
if (debug) std::cout << "End Constructor of AR1" << std::endl;
}

// wrt. parameter_K (bounded parameter)
SparseMatrix<double> AR::getK(const VectorXd& theta_K) const {
    if (type==Type::rw) return C + G;

    assert (theta_K.size() == 1);
    double theta;
    if (type==Type::ar)
        theta = th2a(theta_K(0));

// std::cout << "theta in get K = " << theta << std::endl;
    SparseMatrix<double> K = theta * C + G;
    return K;
}

SparseMatrix<double> AR::get_dK(int index, const VectorXd& alpha) const {
    assert(index==0);
    return C;
}

// compute numerical dK
void AR::update_num_dK() {
    double alpha = th2a(theta_K(0) + eps);
    SparseMatrix<double> K_add_eps = alpha * C + G;
    dK = (K_add_eps - K) / eps;
}

// return length 1 vectorxd : grad_kappa * dkappa/dtheta
VectorXd AR::grad_theta_K() {
    if (numer_grad) return numerical_grad();

    SparseMatrix<double> dK = get_dK_by_index(0);

    double a = th2a(theta_K(0));
    double th = a2th(a);
    double da  = 2 * (exp(th) / pow(1+exp(th), 2));
    // double d2a = 2 * (exp(th) * (-1+exp(th)) / pow(1+exp(th), 3));

    double ret = 0;
    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        VectorXd V = vars[i].getV();
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

        // analytical gradient and numerical hessian
        double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
        double grad = trace - tmp;
        ret += - grad * da / W_size;
    }

    return VectorXd::Constant(1, ret);
    // if (debug) std::cout << "tmp =" << tmp << std::endl;
    // if (debug) std::cout << "trace =" << trace << std::endl;
        // result : trace ~= 0
        //          tmp ~= 20-70
}

void AR::update_each_iter() {
    K = getK(theta_K);
    dK = get_dK(0, theta_K);
    d2K = 0 * C;

    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);

    if (use_num_dK) {
        update_num_dK();
    }

    // trace == 0 for ar case
    // if (!numer_grad && (W_size == V_size))
        // compute_trace();
}

// return length 1 vectorxd : grad_kappa * dkappa/dtheta
VectorXd AR::grad_theta_K(
    SparseMatrix<double>& K,
    SparseMatrix<double>& dK,
    vector<VectorXd>& Ws,
    vector<VectorXd>& prevWs,
    vector<Var>& vars,
    const VectorXd& mu,
    const VectorXd& sigma,
    const VectorXd& h,
    double trace,
    int W_size
) {
    // if (numer_grad) return numerical_grad();

    double a = th2a(theta_K(0));
    double th = a2th(a);
    double da  = 2 * (exp(th) / pow(1+exp(th), 2));
    double ret = 0;
    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        VectorXd V = vars[i].getV();
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
        double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
        double grad = trace - tmp;
        ret += - grad * da / W_size;
    }

    return VectorXd::Constant(1, ret);
}