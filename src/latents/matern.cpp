/*
    Matern model with stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = kappa
        K = kappa^2 * C + G
*/

#include "../latent.h"

Matern::Matern(Rcpp::List& model_list, unsigned long seed)
: Latent(model_list, seed),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
    alpha       (Rcpp::as<int> (model_list["alpha"])),
    Cdiag       (C.diagonal())
{
std::cout << "begin Constructor of Matern " << std::endl;
    symmetricK = true;

    // Init K and Q
    K = getK(parameter_K);
    SparseMatrix<double> Q = K.transpose() * K;

    if (!use_iter_solver) {
        chol_solver_K.init(W_size, 0,0,0);
        chol_solver_K.analyze(K);
    } else {
        CG_solver_K.init(W_size, W_size, W_size, 0.5);
        CG_solver_K.analyze(K);
    }
    compute_trace();
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);

std::cout << "finish Constructor of Matern " << std::endl;
}

SparseMatrix<double> Matern::getK(const VectorXd& parameter_K) const {
    double kappa = parameter_K(0);
    int W_size = G.rows();

    SparseMatrix<double> K_a (W_size, W_size);
        // VectorXd k2C = (kappa * kappa * Cdiag);
        // SparseMatrix<double> KCK = k2C.asDiagonal();
    SparseMatrix<double> KCK = kappa * kappa * C;

    // VectorXd kappas = VectorXd::Constant(W_size, parameter_K(0));
    // SparseMatrix<double> KCK (W_size, W_size);
    //     KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();

    if (alpha==2) {
        // K_a = T (G + KCK) C^(-1/2) -> Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
        K_a = (G + KCK);
    } else if (alpha==4) {
        // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2) -> Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
        K_a = (G + KCK) *
            Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
    } else {
        throw("alpha not equal to 2 or 4 is not implemented");
    }

    return K_a;
}

// stationary
SparseMatrix<double> Matern::get_dK(int index, VectorXd parameter_K) const {
    assert(index==0);
    double kappa = parameter_K(0);
    int W_size = G.rows();
    SparseMatrix<double> dK (W_size, W_size);

    if (alpha==2)
        dK = 2*kappa*C;
    else if (alpha==4)
        dK = 4*kappa*C * G + 4* pow(kappa, 3) * C;
    else
        throw("alpha != 2 or 4");
    return dK;
}

// compute numerical dK
void Matern::update_num_dK() {
    double kappa = parameter_K(0);
    double eps = 0.01;
    SparseMatrix<double> K_add_eps = pow(kappa + eps, 2) * C + G;
    dK = (K_add_eps - K) / eps;
}

VectorXd Matern::get_unbound_theta_K() const {
    assert (parameter_K.size() == 1);

    double th = k2th(parameter_K(0));
    return VectorXd::Constant(1, th);
}

// return length 1 vectorxd : grad_kappa * dkappa/dtheta
VectorXd Matern::grad_theta_K() {
    SparseMatrix<double> dK = get_dK_by_index(0);
    VectorXd V = getV();
    VectorXd SV = getSV();

    VectorXd kappa = parameter_K;
    double th = k2th(kappa(0));

    double da  = exp(th);
    double d2a = exp(th);

    double ret = 0;
    if (numer_grad) {
        // 1. numerical gradient
        ret = numerical_grad()(0);
    } else {
        // 2. analytical gradient and numerical hessian
        double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
        double grad = trace - tmp;

    // sth wrong with hessian?
    // if (debug) std::cout << "tmp =" << tmp << std::endl;
    // if (debug) std::cout << "trace =" << trace << std::endl;

        if (!use_precond) {
            ret = - grad * da / W_size;
        } else {
            // compute numerical hessian
            SparseMatrix<double> K2 = getK_by_eps(0, eps);
            SparseMatrix<double> dK2 = get_dK_by_eps(0, 0, eps);

            // grad(x+eps) - grad(x) / eps
            VectorXd prevV = getPrevV();
            VectorXd prevSV = getPrevSV();
            double grad2_eps = trace_eps - (dK2*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K2 * prevW +  (h - prevV).cwiseProduct(mu));
            double grad_eps  = trace - (dK*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K * prevW +  (h - prevV).cwiseProduct(mu));

            double hess = (grad2_eps - grad_eps) / eps;
    // if (debug) std::cout << "hess =" << hess << std::endl;

            // ret = (grad * da) / (hess * da * da + grad_eps * d2a); reduced to
            ret = grad / (hess * da + grad_eps);
        }
    }

    // if (debug) std::cout << "ret =" << ret << std::endl;
    return VectorXd::Constant(1, ret);
}

void Matern::set_unbound_theta_K(VectorXd theta) {
    double kappa = th2k(theta(0));

    // update theta_K, K and dK
    parameter_K = VectorXd::Constant(1, kappa);
    K = getK(parameter_K);
    dK = get_dK(0, parameter_K);

    if (use_num_dK) {
        update_num_dK();
    }

    if (!numer_grad) compute_trace();
}


// class matern_ope : public Operator {
// private:
//     int alpha;
//     SparseMatrix<double, 0, int> G, C;
//     VectorXd Cdiag;

// public:
//     matern_ope(Rcpp::List& model_list)
//     :   Operator    (model_list),
//         alpha       ( Rcpp::as<int> (model_list["alpha"])),
//         G           ( Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"]) ),
//         C           ( Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"]) ),
//         Cdiag       ( C.diagonal() )
//     {}

//     // set kappa
//     void set_parameter(VectorXd kappa) {
//         assert (kappa.size() == 1);
//         this->parameter_K = kappa;

//         K = getK(kappa);
//         dK = get_dK(0, kappa);

//         if (use_num_dK) {
//             update_num_dK();
//         }
//     }

//     // stationary
//     SparseMatrix<double> get_dK(int index, VectorXd parameter_K) const {
//         assert(index==0);
//         double kappa = parameter_K(0);
//         int W_size = G.rows();
//         SparseMatrix<double> dK (W_size, W_size);

//         if (alpha==2)
//             dK = 2*kappa*C;
//         else if (alpha==4)
//             dK = 4*kappa*C * G + 4* pow(kappa, 3) * C;
//         else
//             throw("alpha != 2 or 4");
//         return dK;
//     }

//     // compute numerical dK
//     void update_num_dK() {
//         double kappa = parameter_K(0);
//         double eps = 0.01;
//         SparseMatrix<double> K_add_eps = pow(kappa + eps, 2) * C + G;
//         dK = (K_add_eps - K) / eps;
//     }
// };