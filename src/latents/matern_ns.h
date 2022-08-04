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
class nonstationaryGC : public Operator {
private:
    int alpha; 
    SparseMatrix<double, 0, int> G, C;
    MatrixXd Bkappa;
    VectorXd Cdiag;
public:
    nonstationaryGC(Rcpp::List ope_in) 
    :   Operator    ( ope_in),
        alpha       ( Rcpp::as<int> (ope_in["alpha"])),
        G           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]) ),
        C           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]) ),
        Bkappa      ( Rcpp::as<MatrixXd> (ope_in["B.kappa"]) ),
        Cdiag       ( C.diagonal() )
    {}

    // here C is diagonal
    void set_parameter(VectorXd theta_kappa) {
        this->parameter_K = theta_kappa;
        
        // update K
        K = getK(theta_kappa);
    }

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
std::cout << "reach here dk1 " << std::endl;
        VectorXd kappas = (Bkappa * parameter_K).array().exp();
std::cout << "reach here dk2 " << std::endl;
        
        int n_dim = G.rows();
        SparseMatrix<double> dK_a (n_dim, n_dim);

        // dKCK
        SparseMatrix<double> CK(n_dim, n_dim);
        CK = kappas.cwiseProduct(Cdiag).asDiagonal();
std::cout << "reach here dk3 " << std::endl;

        SparseMatrix<double> dKCK(n_dim, n_dim);
        dKCK = kappas.cwiseProduct(Bkappa.col(index)).asDiagonal() * CK + CK * kappas.cwiseProduct(Bkappa.col(index)).asDiagonal(); // kappas * (Bkappa * CK + CK * Bkappa).sparseView();
std::cout << "reach here dk4 " << std::endl;

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
std::cout << "reach here dk5 " << std::endl;

        return dK_a;
    }

};

// get_K_params, grad_K_params, set_K_params, output
class Matern_ns : public Latent {
public:
    Matern_ns(Rcpp::List latent_in) 
    : Latent(latent_in)
    {
if (debug) std::cout << "constructor of matern ns" << std::endl;
        symmetricK = true;
        Rcpp::List operator_in = Rcpp::as<Rcpp::List> (latent_in["operator_in"]); // containing C and G
        ope = new nonstationaryGC(operator_in);
            Rcpp::List start = Rcpp::as<Rcpp::List> (latent_in["start"]);
            VectorXd parameter_K = Rcpp::as< VectorXd > (start["theta_K"]);
            ope->set_parameter(parameter_K);
        
        // Init K and Q
        SparseMatrix<double> K = getK();
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

    // generating output
    Rcpp::List get_estimates() const {
        return Rcpp::List::create(
            Rcpp::Named("theta.kappa") = ope->get_parameter(),
            Rcpp::Named("theta.mu")    = theta_mu,
            Rcpp::Named("theta.sigma") = theta_sigma,
            Rcpp::Named("theta.noise") = var->get_var()
        );
    }

    VectorXd grad_theta_K() {
std::cout << "begin grad theta K " << std::endl;
        n_ope = ope->get_n_params();
std::cout << "n_ope = " << n_ope << std::endl;

        SparseMatrix<double> K = ope->getK();
        // SparseMatrix<double> dK = ope->get_dK(0);        
        // SparseMatrix<double> dK (n_mesh, n_mesh);
        // for (int i=0; i < n_ope; i++) {
        //     dK = dK + ope->get_dK(i);
        // }

        VectorXd V = getV();
        VectorXd SV = getSV();
        
        VectorXd kappa = ope->get_parameter();

        VectorXd grad;
        if (numer_grad) {
            // 1. numerical gradient
            grad = numerical_grad();
        } else { 
            // 2. analytical gradient and numerical hessian
            for (int i=0; i < n_ope; i++) {
                // dK for each index
                SparseMatrix<double> dK = ope->get_dK(i);
                
                VectorXd dkw = dK*W;
std::cout << "size= " << dkw.size() << std::endl;

                VectorXd tmp2 = K * W + (h - V).cwiseProduct(mu);
std::cout << "reach here 1 " << std::endl;
                double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(tmp2);
std::cout << "reach here 2 " << std::endl;

                // compute trace
                if (!symmetricK) {
                    lu_solver_K.computeKTK(K);
                    trace = lu_solver_K.trace(dK);
                } else {
                    chol_solver_K.compute(K);
                    trace = chol_solver_K.trace(dK);
                }

                grad(i) = trace - tmp;
            }
        }

std::cout << "grad= " << grad << std::endl;
        return grad;
    }

};

#endif