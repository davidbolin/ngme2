#include "../operator.h"

General::General(const Rcpp::List& operator_list):
    Operator(operator_list),
    theta_trans(Rcpp::as<std::string>(operator_list["theta_trans"])),
    matrices(Rcpp::as<std::vector<SparseMatrix<double, 0, int>>>(operator_list["matrices"]))
{}

void General::update_K(const VectorXd& theta_K) {
    // update K using theta_K
    VectorXd theta = param_trans_fun(theta_K, theta_trans);
    K.setZero();
    for (int i = 0; i < matrices.size(); i++) {
        K += theta[i] * matrices[i];
    }
}

void General::update_dK(const VectorXd& theta_K) {
    // update dK using theta_K
}

VectorXd General::param_trans_fun(const VectorXd& theta_K, const std::string& name) const {
    // transform theta_K using trans_name
    if (name == "matern") {
        Vector2d trans_theta = Vector2d::Zero();    
        trans_theta(0) = exp(2 * theta_K(0));
        trans_theta(1) = 1;
        return trans_theta;
    } else {
        return theta_K;
    }
}
