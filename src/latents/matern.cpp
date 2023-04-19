/*
    Matern model with stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = kappa
        K = kappa^2 * C + G
*/

#include "../latent.h"

Matern::Matern(const Rcpp::List& model_list, unsigned long seed)
: Latent(model_list, seed),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
    alpha       (Rcpp::as<int> (model_list["alpha"])),
    Cdiag       (C.diagonal())
{
// std::cout << "begin Constructor of Matern " << std::endl;
    symmetricK = true;
    // Init K and Q
    K = getK(theta_K);
    chol_solver_K.init(W_size, 0,0,0);
    chol_solver_K.analyze(K);

    SparseMatrix<double> Q = K.transpose() * K;
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);

    update_each_iter();
std::cout << "finish Constructor of Matern " << std::endl;
}

SparseMatrix<double> Matern::getK(const VectorXd& theta_K) const {
    double kappa = exp(theta_K(0));
    int W_size = G.rows();

    SparseMatrix<double> K_a (W_size, W_size);
        // VectorXd k2C = (kappa * kappa * Cdiag);
        // SparseMatrix<double> KCK = k2C.asDiagonal();
    SparseMatrix<double> KCK = kappa * kappa * C;

    // VectorXd kappas = VectorXd::Constant(W_size, theta_K(0));
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
SparseMatrix<double> Matern::get_dK(int index, const VectorXd& theta_K) const {
    assert(index==0);
    double kappa = exp(theta_K(0));
    int W_size = G.rows();
    SparseMatrix<double> dK (W_size, W_size);

    if (alpha==2)
        dK = 2.0*kappa*C;
    else if (alpha==4)
        dK = 4.0*kappa*C * G + 4.0* pow(kappa, 3) * C;
    else
        throw("alpha != 2 or 4");

    // dkappa / dtheta = kappa
    return kappa * dK;
}

// ---------------------------- Matern_ns ----------------------------

/*
    Matern model with non-stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = theta.kappa
        kappas =  exp(Bkappa * theta.kappa)
        K = kappa^2 * C + G
*/
Matern_ns::Matern_ns(const Rcpp::List& model_list, unsigned long seed, Type type)
: Latent(model_list, seed),
    type        (type),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
    alpha       (2),
    Bkappa      (Rcpp::as<MatrixXd> (model_list["B_K"])),
    Cdiag       (C.diagonal())
{
//  std::cout << "constructor of matern ns" << std::endl;
    K = getK(theta_K);

    if (type==Type::matern_ns) {
        alpha = Rcpp::as<int> (model_list["alpha"]);
        symmetricK = true;
        chol_solver_K.init(W_size, 0,0,0);
        chol_solver_K.analyze(K);
    } else {
        symmetricK = false;
        lu_solver_K.init(W_size, 0,0,0);
        lu_solver_K.analyze(K);
    }

    // Init Q
    SparseMatrix<double> Q = K.transpose() * K;
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);

    update_each_iter();
//  std::cout << "finish constructor of matern ns" << std::endl;
}

// inherit get_K_parameter, grad_K_parameter, set_K_parameter

SparseMatrix<double> Matern_ns::getK(const VectorXd& theta_kappa) const {
    VectorXd kappas = (Bkappa * theta_kappa).array().exp();
    // std::cout <<  "theta_kappa here = " << theta_kappa << std::endl;

    int n_dim = G.rows();
    SparseMatrix<double> K_a (n_dim, n_dim);
    if (type == Type::matern_ns) {
        SparseMatrix<double> KCK (n_dim, n_dim);
            KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
        if (alpha==2) {
            // K_a = T (G + KCK) C^(-1/2)
            // Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
            K_a = (G + KCK);
        } else if (alpha==4) {
            // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2)
            // Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
            K_a = (G + KCK) * Cdiag.cwiseInverse().asDiagonal() *
            (G + KCK);
        } else {
            throw("alpha not equal to 2 or 4 is not implemented");
        }
    } else if (type == Type::ou) {
        K_a = kappas.asDiagonal() * C + G;
    }

    return K_a;
}

// dK wrt. theta_K[index]
SparseMatrix<double> Matern_ns::get_dK(int index, const VectorXd& params) const {
    VectorXd kappas = (Bkappa * theta_K).array().exp();

    int n_dim = G.rows();
    SparseMatrix<double> dK_a (n_dim, n_dim);

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
            dK_a = dKCK;
        }
        else if (alpha == 4) {
            SparseMatrix<double> KCK(n_dim, n_dim);
            KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
            SparseMatrix<double> tmp = Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
            dK_a = dKCK * tmp + tmp * dKCK;
        }
        else {
            throw("alpha not equal to 2 or 4 is not implemented");
        }
    } else {
        // check
        dK_a = kappas.cwiseProduct(Bkappa.col(index)).asDiagonal() * C;
    }

    return dK_a;
}