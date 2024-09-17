#include "../operator.h"

General::General(const Rcpp::List& operator_list):
    Operator(operator_list),
    matrices(Rcpp::as<std::vector<SparseMatrix<double, 0, int>>>(operator_list["matrices"])),
    trans_type(Rcpp::as<std::string>(operator_list["trans_type"]))
{}

void General::update_K(const VectorXd& theta_K) {
    VectorXd theta = param_trans_fun(theta_K, trans_type);
    K.setZero();
    for (int i = 0; i < matrices.size(); i++) {
        K += theta[i] * matrices[i];
    }
}

void General::update_dK(const VectorXd& theta_K) {
    // update dK using theta_K
}

VectorXd General::param_trans_fun(const VectorXd& theta_K, const std::string& trans_type) const {
    //  -1 + (2 * exp(th)) / (1 + exp(th))
    auto th2a = [](double th) {
        return -1 + (2 * exp(th)) / (1 + exp(th));
    };

    // transform theta_K using trans_name
    if (trans_type == "matern") {
        Vector2d trans_theta = Vector2d::Zero();    
        trans_theta(0) = exp(2 * theta_K(0));
        trans_theta(1) = 1;
        return trans_theta;
    } else if (trans_type == "ar1 x matern") {
        Vector4d trans_theta = Vector4d::Zero();
        trans_theta(0) = th2a(theta_K(0)) * exp(2 * theta_K(1));
        trans_theta(1) = th2a(theta_K(0));
        trans_theta(2) = exp(2 * theta_K(1));
        trans_theta(3) = 1;
        return trans_theta;
    } else {
        return theta_K;
    }
}
