/*
    Matern model with stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = kappa
        K = kappa^2 * C + G
*/

#include "../operator.h"

Matern::Matern(const Rcpp::List& operator_list):
    Operator(operator_list),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (operator_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (operator_list["C"])),
    alpha       (Rcpp::as<int> (operator_list["alpha"])),
    Cdiag       (C.diagonal())
{}

void Matern::update_K(const VectorXd& theta_K) {
    double kappa = exp(theta_K(0));
    int W_size = G.rows();
    SparseMatrix<double> KCK = kappa * kappa * C;

    // for test old ngme model parameterization
    // if (theta_K.size() == 1 && alpha == 2) {
    //     // stationary
    //     K = 0.5 * (pow(kappa, -1.5) * G + pow(kappa, 0.5) * C);
    // } else
    if (alpha==2) {
        // K_a = T (G + KCK) C^(-1/2) -> Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
        K = (G + KCK);
    } else if (alpha==4) {
        // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2) -> Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
        K = (G + KCK) *
            Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
    } else {
        throw("alpha not equal to 2 or 4 is not implemented");
    }
// std::cout << " K = " << K << std::endl;
}

void Matern::update_dK(const VectorXd& theta_K) {
    assert(theta_K.size()==1);
    double kappa = exp(theta_K(0));
    int W_size = G.rows();

    if (alpha==2)
        dK[0] = 2.0*kappa*C;
    else if (alpha==4)
        dK[0] = 4.0*kappa*C * G + 4.0* pow(kappa, 3) * C;
    else
        throw("alpha != 2 or 4");

    // dkappa / dtheta = kappa
    dK[0] = kappa * dK[0];
}

// ---------------------------- Matern_ns ----------------------------

/*
    Matern model with non-stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = theta.kappa
        kappas =  exp(Bkappa * theta.kappa)
        K = kappa^2 * C + G
*/
Matern_ns::Matern_ns(const Rcpp::List& operator_list, Type type):
    Operator(operator_list),
    type        (type),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (operator_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (operator_list["C"])),
    alpha       (2),
    Bkappa      (Rcpp::as<MatrixXd> (operator_list["B_K"])),
    Cdiag       (C.diagonal())
{
//  std::cout << "constructor of matern ns" << std::endl;
}

// inherit get_K_parameter, grad_K_parameter, set_K_parameter
void Matern_ns::update_K(const VectorXd& theta_kappa) {
    VectorXd kappas = (Bkappa * theta_kappa).array().exp();
    // std::cout <<  "theta_kappa here = " << theta_kappa << std::endl;

    int n_dim = G.rows();
    if (type == Type::matern_ns) {
        SparseMatrix<double> KCK (n_dim, n_dim);
            KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
        if (alpha==2) {
            // K_a = T (G + KCK) C^(-1/2)
            // Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
            K = (G + KCK);
        } else if (alpha==4) {
            // K = T (G + KCK) C^(-1) (G+KCK) C^(-1/2)
            // Actually, K = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
            K = (G + KCK) * Cdiag.cwiseInverse().asDiagonal() *
            (G + KCK);
        } else {
            throw("alpha not equal to 2 or 4 is not implemented");
        }
    } else if (type == Type::ou) {
        K = kappas.asDiagonal() * C + G;
    }
}

// dK wrt. theta_K[index]
void Matern_ns::update_dK(const VectorXd& theta_K) {
// std::cout << "update_dK matern" << std::endl;
    VectorXd kappas = (Bkappa * theta_K).array().exp();
    int n_dim = G.rows();

    for (int index = 0; index < n_theta_K; index++) {
        if (type == Type::matern_ns) {
            // dKCK
            SparseMatrix<double> CK(n_dim, n_dim);
            // CK = kappas.cwiseProduct(Cdiag).asDiagonal();

            SparseMatrix<double> dKCK(n_dim, n_dim);
            //  dKCK = 2*kappas.cwiseProduct(Bkappa.col(index)).asDiagonal() * CK;
            // kappas * (Bkappa * CK + CK * Bkappa).sparseView();
                VectorXd kappas2 = kappas.cwiseProduct(kappas);
                dKCK = 2*kappas2.cwiseProduct(Cdiag).cwiseProduct(Bkappa.col(index)).asDiagonal();
            if (alpha == 2) {
                dK[index] = dKCK;
            }
            else if (alpha == 4) {
                SparseMatrix<double> KCK(n_dim, n_dim);
                KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
                SparseMatrix<double> tmp = Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
                dK[index] = dKCK * tmp + tmp * dKCK;
            }
            else {
                throw("alpha not equal to 2 or 4 is not implemented");
            }
        } else {
            // check
            dK[index] = kappas.cwiseProduct(Bkappa.col(index)).asDiagonal() * C;
        }
    }
}