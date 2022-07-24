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

    SparseMatrix<double> get_dK(int index, VectorXd params) const {
        return G;        

        // std::vector<SparseMatrix<double>> dK_a(Bkappa.cols());
        // for (int i = 0; i < Bkappa.cols(); i++)
        // {
        //     VectorXd kappas = (Bkappa * parameter_K).array().exp();

        //     int n_dim = G.rows();

        //     // dKCK
        //     SparseMatrix<double> CK(n_dim, n_dim);
        //     CK = kappas.cwiseProduct(Cdiag).asDiagonal();

        //     SparseMatrix<double> dKCK(n_dim, n_dim);
        //     dKCK = kappas.cwiseProduct(Bkappa.col(i)) * CK + CK * kappas.cwiseProduct(Bkappa.col(i)); // kappas * (Bkappa * CK + CK * Bkappa).sparseView();

        //     if (alpha == 2)
        //     {
        //         dK_a[i] = dKCK;
        //     }
        //     else if (alpha == 4)
        //     {
        //         SparseMatrix<double> KCK(n_dim, n_dim);
        //         KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
        //         SparseMatrix<double> tmp = Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
        //         dK_a[i] = dKCK * tmp + tmp * dKCK;
        //     }
        //     else
        //     {
        //         throw("alpha not equal to 2 or 4 is not implemented");
        //     }
        // }

        // return dK_a;
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
};

#endif