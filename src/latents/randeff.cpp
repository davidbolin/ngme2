// Notice here:
// for R interface, theta_K is directly the parameter_K of the Opeartor object
// for C interface, theta_K is f(parameter_K) such that it is unbounded
#include "../operator.h"

/*
    Random effect model:
        parameter_K = invSigma
*/

// W_size = V_size
// get_K_params, grad_K_params, set_K_params, output
Randeff::Randeff(const Rcpp::List& operator_list):
    Operator(operator_list),
    n_repl(operator_list["n_repl"])
{
std::cout << "End Constructor of Randeff1" << std::endl;
}

//     VectorXd getM() const {
//         double V = var.getV()(0);
// V = 1;
//         return V*(V-1) * Sigma * mu;
//     }

// K^T K = Sigma^-1
void Randeff::update_K(const VectorXd& theta_K) {
    VectorXd diag = theta_K.head(W_size).array().exp();
    VectorXd offdiag = theta_K.tail(W_size * (W_size+1) / 2 - W_size);

    int index = 0;
    MatrixXd L (W_size, W_size); L.setZero();
    L.diagonal() = diag;
    for (int col=0; col < W_size; col++) {
        for (int row=col+1; row < W_size; row++) {
            L(row, col) = offdiag(index);
            index++;
        }
    }

    K = L.sparseView();
}

MatrixXd Randeff::get_dK_dense(int index, const VectorXd& alpha) const {
    // MatrixXd dK = MatrixXd::Zero(W_size, W_size);
    // VectorXd tmp = VectorXd::Zero(W_size * (W_size+1) / 2);
    // if (index < W_size) {
    //     dK(index, index) = exp(alpha(index));
    // } else {
    //     tmp(index) = 1;
    //     dK = getK(tmp);
    //     dK.diagonal().setZero();
    // }
    // return dK;
}

void Randeff::update_dK(const VectorXd& theta_K) {

}

// SparseMatrix<double> Randeff::get_dK(int index, const VectorXd& alpha) const {
//     return get_dK_dense(index, alpha).sparseView();
// }

// // handle specially, always return 0, but update Sigma
// VectorXd Randeff::grad_theta_K() {
//     VectorXd grad = VectorXd::Zero(n_theta_K);
//     MatrixXd K = getK(theta_K).toDense();

//     for (int i=0; i < n_rep; i++) {
//         double V = vars[i].getV()(0);
//         VectorXd W = Ws[i];
//         VectorXd tmp = 1/V * (K*W + (1-V)*mu);
//         for (int j=0; j < n_theta_K; j++) {
//             MatrixXd dK = get_dK_dense(j, theta_K);
// // std::cout << "dK \n" << dK << std::endl;
//             if (j < W_size) {
//                 grad(j) = K.llt().solve(dK).diagonal().sum() -
//                     W.transpose() * dK.transpose() * tmp;
//                 grad(j) = 1.0/sqrt(n_repl) * grad(j);
//             } else {
//                 // how to choose the step size?
//                 grad(j) = -1.0/sqrt(n_repl) * W.transpose() * dK.transpose() * tmp;
//             }
//         }
//     }
//     grad = - grad / W_size;
// // std::cout << "grad_theta_K = " << grad << std::endl;
//     return grad;
// }


// void Randeff::sample_cond_V() {
// // return Operator::sample_cond_V();
//     VectorXd a_inc_vec = mu.cwiseQuotient(sigma).array().pow(2);
//     for (int i=0; i < n_rep; i++) {
//         VectorXd W = Ws[i];
//         VectorXd tmp = (K * W + mu.cwiseProduct(h));
//         VectorXd b_inc_vec = tmp.cwiseQuotient(sigma).array().pow(2);
//         vars[i].sample_cond_V(a_inc_vec, b_inc_vec, W_size, true);
//     }
// }

// const VectorXd Randeff::get_grad() {
//     VectorXd grad = Operator::get_grad();
// // std::cout << "V = " << vars[0].getV() << std::endl;
// // std::cout << "mu = " << mu << std::endl;
// // std::cout << "grad = " << grad << std::endl;
//     return grad;
// }

// void Randeff::update_each_iter() {
//         mu = (B_mu * theta_mu);
//         sigma = (B_sigma * theta_sigma).array().exp();
//         if (noise_type=="normal_nig") sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

//         // update on K, dK, trace, bigK for sampling...
//         K = getK(theta_K);

//         for (int i=0; i < n_rep; i++)
//             setSparseBlock(&K_rep, i*V_size, i*V_size, K);
//     }