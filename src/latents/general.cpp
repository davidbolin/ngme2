#include "../operator.h"

General::General(const Rcpp::List& operator_list):
    Operator(operator_list),
    matrices(Rcpp::as<std::vector<SparseMatrix<double, 0, int>>>(operator_list["matrices"])),
    idx_mat(Rcpp::as<MatrixXd>(operator_list["idx_mat"])),
    trans(Rcpp::as<std::vector<std::string>>(operator_list["trans"]))
{}

void General::update_K(const VectorXd& theta_K) {
    VectorXd coef = compute_coef(theta_K, idx_mat, trans);
    K.setZero();
    for (int i = 0; i < matrices.size(); i++) {
        K += coef[i] * matrices[i];
    }
}

void General::update_dK(const VectorXd& theta_K) {
    // update dK using theta_K
}

// VectorXd General::param_trans_fun(const VectorXd& theta_K, const std::string& trans) const {
//     //  -1 + (2 * exp(th)) / (1 + exp(th))
//     auto th2a = [](double th) {
//         return -1 + (2 * exp(th)) / (1 + exp(th));
//     };

//     // transform theta_K using trans_name
//     if (trans == "matern") {
//         Vector2d trans_theta = Vector2d::Zero();    
//         trans_theta(0) = exp(2 * theta_K(0));
//         trans_theta(1) = 1;
//         return trans_theta;
//     } else if (trans_type == "ar1 x matern") {
//         Vector4d trans_theta = Vector4d::Zero();
//         trans_theta(0) = th2a(theta_K(0)) * exp(2 * theta_K(1));
//         trans_theta(1) = th2a(theta_K(0));
//         trans_theta(2) = exp(2 * theta_K(1));
//         trans_theta(3) = 1;
//         return trans_theta;
//     } else {
//         return theta_K;
//     }
// }


// template of idx_mat
//        rho  beta kappa
// [1,]  TRUE FALSE  TRUE
// [2,] FALSE  TRUE  TRUE
// [3,]  TRUE FALSE FALSE
// [4,] FALSE  TRUE FALSE


VectorXd General::compute_coef(const VectorXd& theta_K, const MatrixXd& idx_mat, const std::vector<std::string>& trans) const {
    int nrow = idx_mat.rows();
    VectorXd coef = VectorXd::Ones(nrow);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < idx_mat.cols(); j++) {
            if (idx_mat(i, j) == 1) {
                if (trans[j] == "exp2") {  
                    coef[i] *= exp(2 * theta_K[j]);
                } else if (trans[j] == "tanh") {
                    double th = theta_K[j];
                    coef[i] *= (-1 + (2 * exp(th)) / (1 + exp(th)));
                } else if (trans[j] == "identity") {
                    coef[i] *= theta_K[j];
                }
            }
        }
    }
    return coef;
}