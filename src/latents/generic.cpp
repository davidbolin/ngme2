#include "../operator.h"

Generic::Generic(const Rcpp::List& operator_list):
    Operator(operator_list),
    matrices(Rcpp::as<std::vector<SparseMatrix<double, 0, int>>>(operator_list["matrices"])),
    idx_mat(Rcpp::as<MatrixXd>(operator_list["idx_mat"])),
    trans(Rcpp::as<std::vector<std::string>>(operator_list["trans"]))
{}

void Generic::update_K(const VectorXd& theta_K) {
    VectorXd coef = compute_coef(theta_K, idx_mat, trans);
    K.setZero();
    for (int i = 0; i < matrices.size(); i++) {
        K += coef[i] * matrices[i];
    }
}

void Generic::update_dK(const VectorXd& theta_K) {
    // update dK using theta_K
}

VectorXd Generic::compute_coef(const VectorXd& theta_K, const MatrixXd& idx_mat, const std::vector<std::string>& trans) const {
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