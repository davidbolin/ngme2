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
    n_reff(K.rows())
{
// std::cout << "End Constructor of Randeff1" << std::endl;
}

//     VectorXd getM() const {
//         double V = var.getV()(0);
// V = 1;
//         return V*(V-1) * Sigma * mu;
//     }

// K^T K = Sigma^-1
void Randeff::update_K(const VectorXd& theta_K) {
    VectorXd diag = theta_K.head(n_reff).array().exp();
    VectorXd offdiag = theta_K.tail(n_reff * (n_reff+1) / 2 - n_reff);

    int index = 0;
    MatrixXd L (n_reff, n_reff); L.setZero();
    L.diagonal() = diag;
    for (int col=0; col < n_reff; col++) {
        for (int row=col+1; row < n_reff; row++) {
            L(row, col) = offdiag(index);
            index++;
        }
    }

    K = L.sparseView();
    // std::cout << "K \n" << K << std::endl;
}

void Randeff::update_dK(const VectorXd& theta_K) {
    // std::cout << "Start update_dK" << std::endl;
    for (int index=0; index < n_theta_K; index++) {
        if (index < n_reff) {
            dK[index].coeffRef(index, index) = exp(theta_K(index));
        } else {
            // compute row and col given index
            int idx = index - n_reff;
            int n = n_reff - 1;
            int tmp = n;
            int col = 0;
            while (idx >= tmp) {
                n--; col++;
                tmp += n;
            }
            int row = n_reff - (tmp - idx);
            dK[index].coeffRef(row, col) = 1;
        }
    // std::cout << "dK[" << index << "] \n" << dK[index] << std::endl;
    }
    // std::cout << "end update_dK" << std::endl;
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