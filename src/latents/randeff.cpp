// Notice here:
// for R interface, theta_K is directly the parameter_K of the Opeartor object
// for C interface, theta_K is f(parameter_K) such that it is unbounded
#include "../latent.h"

/*
    Random effect model:
        parameter_K = invSigma
*/

// W_size = V_size
// get_K_params, grad_K_params, set_K_params, output
Randeff::Randeff(const Rcpp::List& model_list, unsigned long seed)
    : Latent(model_list, seed)
{
if (debug) std::cout << "Begin Constructor of Randeff1" << std::endl;

    // Init K and Q
    K = getK(theta_K);
// if (debug) std::cout << "Begin Constructor of Randeff2" << std::endl;
    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);
// if (debug) std::cout << "K = " << K << std::endl;

    // Init Q
    SparseMatrix<double> Q = K.transpose() * K;
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);
if (debug) std::cout << "End Constructor of Randeff1" << std::endl;
}

//     VectorXd getM() const {
//         double V = var.getV()(0);
// V = 1;
//         return V*(V-1) * Sigma * mu;
//     }

// K^T K = Sigma^-1
SparseMatrix<double> Randeff::getK(const VectorXd& theta_K) const {
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

    SparseMatrix<double> K = L.sparseView();
    return K;
}

MatrixXd Randeff::get_dK_dense(int index, const VectorXd& alpha) const {
    MatrixXd dK = MatrixXd::Zero(W_size, W_size);
    VectorXd tmp = VectorXd::Zero(W_size * (W_size+1) / 2);
    if (index < W_size) {
        dK(index, index) = exp(alpha(index));
    } else {
        tmp(index) = 1;
        dK = getK(tmp);
        dK.diagonal().setZero();
    }
    return dK;
}

SparseMatrix<double> Randeff::get_dK(int index, const VectorXd& alpha) const {
    return get_dK_dense(index, alpha).sparseView();
}

// handle specially, always return 0, but update Sigma
VectorXd Randeff::grad_theta_K() {
    VectorXd grad = VectorXd::Zero(n_theta_K);
    MatrixXd K = getK(theta_K).toDense();
    if (numer_grad) {
        // compute func_K(K)
        double val=0;
        for (int i = 0; i < n_rep; i++) {
            VectorXd W = Ws[i];
            double V = vars[i].getV()(0);
            VectorXd tmp = K * W + (1-V) * mu;
            val += K.diagonal().sum() - 0.5 * (1/V) * tmp.dot(tmp);
        }
        // compute func_K(K_add_eps)
        for (int i=0; i < n_theta_K; i++) {
            MatrixXd K_add_eps = getK_by_eps(i, eps).toDense();
        // std::cout << "K_add_eps \n" << K_add_eps << std::endl;
            double val_add_eps = 0;
            for (int i = 0; i < n_rep; i++) {
                VectorXd W = Ws[i];
                double V = vars[i].getV()(0);
                VectorXd tmp = K_add_eps * W + (1-V) * mu;
                val_add_eps += K_add_eps.diagonal().sum() - 0.5 * (1/V) * tmp.dot(tmp);
            }
            double num_grad = (val_add_eps - val) / eps;
            grad(i) = num_grad / W_size;
        }
    } else {
        // analytical gradient
        for (int i=0; i < n_rep; i++) {
            double V = vars[i].getV()(0);
            VectorXd W = Ws[i];
            VectorXd tmp = 1/V * (K*W + (1-V)*mu);
            for (int j=0; j < n_theta_K; j++) {
                MatrixXd dK = get_dK_dense(j, theta_K);
// std::cout << "dK \n" << dK << std::endl;
                if (j < W_size) {
                    grad(j) = K.llt().solve(dK).diagonal().sum() -
                        W.transpose() * dK.transpose() * tmp;
                } else {
                    // how to choose the step size?
                    grad(j) = -0.05 * W.transpose() * dK.transpose() * tmp;
                }
            }
        }
        grad = - grad / W_size;
    }
// std::cout << "grad_theta_K = " << grad << std::endl;
    return grad;
}

// called after set parameter
void Randeff::update_each_iter() {
    K = getK(theta_K);
    // dK = get_dK(0, theta_K);
    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);
// std::cout << "finish update_each_iter in randeff" << std::endl;
}

void Randeff::sample_cond_V() {
    // std::cout << "sample_cond_V in randeff" << std::endl;
    MatrixXd Sigma = K * K.transpose();
    double a_inc_vec = mu.dot(Sigma.llt().solve(mu));

    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        double b_inc_vec = (W+mu).dot(Sigma.llt().solve(W+mu));
        vars[i].sample_cond_V(a_inc_vec, b_inc_vec, W_size);
    }

    // std::cout << "var.getV = " << vars[0].getV() << std::endl;
}