#ifndef NGME_MATERN_NS_H
#define NGME_MATERN_NS_H

#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;
/*
    Matern model with non-stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = theta.kappa
        kappas =  exp(Bkappa * theta.kappa)
        K = kappa^2 * C + G
*/

// get_K_params, grad_K_params, set_K_params, output
class Matern_ns : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    MatrixXd Bkappa;
    VectorXd Cdiag;
public:
    Matern_ns(Rcpp::List& model_list, unsigned long seed)
    : Latent(model_list, seed),
      G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
      C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
      alpha       (Rcpp::as<int> (model_list["alpha"])),
      Bkappa      (Rcpp::as<MatrixXd> (model_list["B.kappa"])),
      Cdiag       (C.diagonal())
    {
if (debug) std::cout << "constructor of matern ns" << std::endl;
        symmetricK = true;

        // Init K and Q
        K = getK(parameter_K);
        SparseMatrix<double> Q = K.transpose() * K;

        chol_solver_K.init(n_mesh, 0,0,0);
        chol_solver_K.analyze(K);
        // compute_trace();

        // Init Q
        solver_Q.init(n_mesh, 0,0,0);
        solver_Q.analyze(Q);
if (debug) std::cout << "finish constructor of matern ns" << std::endl;
    }

    // inherit get_K_parameter, grad_K_parameter, set_K_parameter

    SparseMatrix<double> getK(VectorXd theta_kappa) const {
        VectorXd kappas = (Bkappa * theta_kappa).array().exp();

        int n_dim = G.rows();
        SparseMatrix<double> K_a (n_dim, n_dim);
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

        return K_a;
    }

    // dK wrt. theta_K[index]
    SparseMatrix<double> get_dK(int index, VectorXd params) const {
        VectorXd kappas = (Bkappa * parameter_K).array().exp();

        int n_dim = G.rows();
        SparseMatrix<double> dK_a (n_dim, n_dim);

        // dKCK
        SparseMatrix<double> CK(n_dim, n_dim);
        // CK = kappas.cwiseProduct(Cdiag).asDiagonal();

        SparseMatrix<double> dKCK(n_dim, n_dim);
        //  dKCK = 2*kappas.cwiseProduct(Bkappa.col(index)).asDiagonal() * CK;
        // kappas * (Bkappa * CK + CK * Bkappa).sparseView();
            VectorXd kappas2 = kappas.cwiseProduct(kappas);
            dKCK = 2*kappas2.cwiseProduct(Cdiag).cwiseProduct(Bkappa.col(index)).asDiagonal();


        if (alpha == 2)
        {
            dK_a = dKCK;
        }
        else if (alpha == 4)
        {
            SparseMatrix<double> KCK(n_dim, n_dim);
            KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
            SparseMatrix<double> tmp = Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
            dK_a = dKCK * tmp + tmp * dKCK;
        }
        else
        {
            throw("alpha not equal to 2 or 4 is not implemented");
        }

        return dK_a;
    }

    VectorXd grad_theta_K() {
std::cout << "matern_ns : begin grad theta K " << std::endl;
        VectorXd V = getV();
        VectorXd SV = getSV();

        VectorXd grad (n_theta_K);
        if (numer_grad) {
            // 1. numerical gradient
            grad = numerical_grad();
        } else {
            // 2. analytical gradient and numerical hessian
            chol_solver_K.compute(K);
            for (int i=0; i < n_theta_K; i++) {
                // dK for each index
                SparseMatrix<double> dK = get_dK_by_index(i);

                VectorXd tmp2 = K * W + (h - V).cwiseProduct(mu);
                double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(tmp2);

                // compute trace
                if (i > 0) {
                    trace = chol_solver_K.trace(dK);
                }

                grad(i) = (trace - tmp) / n_mesh;
            }
        }

std::cout << "matern_ns : grad= " << grad << std::endl;
        return grad;
    }

    void set_unbound_theta_K(VectorXd theta_K) {
std::cout << "matern_ns : begin set theta K " << std::endl;

        parameter_K = theta_K;
        K = getK(parameter_K);
        if (!numer_grad) compute_trace();
    }

    // generating output
    Rcpp::List get_estimates() const {
        return Rcpp::List::create(
            Rcpp::Named("theta.kappa") = parameter_K,
            Rcpp::Named("theta.mu")    = theta_mu,
            Rcpp::Named("theta.sigma") = theta_sigma,
            Rcpp::Named("theta.noise") = var->get_var()
        );
    }

};

#endif

// class nonstationaryGC : public Operator {
// private:
//     int alpha;
//     SparseMatrix<double, 0, int> G, C;
//     MatrixXd Bkappa;
//     VectorXd Cdiag;
// public:
//     nonstationaryGC(Rcpp::List& model_list)
//     :   Operator    (model_list),
//         alpha       ( Rcpp::as<int> (model_list["alpha"])),
//         G           ( Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"]) ),
//         C           ( Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"]) ),
//         Bkappa      ( Rcpp::as<MatrixXd> (model_list["B.kappa"]) ),
//         Cdiag       ( C.diagonal() )
//     {}

//     // here C is diagonal
//     void set_parameter(VectorXd theta_kappa) {
//         this->parameter_K = theta_kappa;

//         // update K
//         K = getK(theta_kappa);
//     }

//     SparseMatrix<double> getK(VectorXd theta_kappa) const {
//         VectorXd kappas = (Bkappa * theta_kappa).array().exp();

//         int n_dim = G.rows();
//         SparseMatrix<double> K_a (n_dim, n_dim);
//         SparseMatrix<double> KCK (n_dim, n_dim);
//             KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();

//         if (alpha==2) {
//             // K_a = T (G + KCK) C^(-1/2)
//             // Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
//             K_a = (G + KCK);
//         } else if (alpha==4) {
//             // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2)
//             // Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
//             K_a = (G + KCK) * Cdiag.cwiseInverse().asDiagonal() *
//             (G + KCK);
//         } else {
//             throw("alpha not equal to 2 or 4 is not implemented");
//         }

//         return K_a;
//     }

//     // dK wrt. theta_K[index]
//     SparseMatrix<double> get_dK(int index, VectorXd params) const {
//         VectorXd kappas = (Bkappa * parameter_K).array().exp();

//         int n_dim = G.rows();
//         SparseMatrix<double> dK_a (n_dim, n_dim);

//         // dKCK
//         SparseMatrix<double> CK(n_dim, n_dim);
//         // CK = kappas.cwiseProduct(Cdiag).asDiagonal();

//         SparseMatrix<double> dKCK(n_dim, n_dim);
//         //  dKCK = 2*kappas.cwiseProduct(Bkappa.col(index)).asDiagonal() * CK;
//         // kappas * (Bkappa * CK + CK * Bkappa).sparseView();
//             VectorXd kappas2 = kappas.cwiseProduct(kappas);
//             dKCK = 2*kappas2.cwiseProduct(Cdiag).cwiseProduct(Bkappa.col(index)).asDiagonal();


//         if (alpha == 2)
//         {
//             dK_a = dKCK;
//         }
//         else if (alpha == 4)
//         {
//             SparseMatrix<double> KCK(n_dim, n_dim);
//             KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
//             SparseMatrix<double> tmp = Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
//             dK_a = dKCK * tmp + tmp * dKCK;
//         }
//         else
//         {
//             throw("alpha not equal to 2 or 4 is not implemented");
//         }

//         return dK_a;
//     }

// };