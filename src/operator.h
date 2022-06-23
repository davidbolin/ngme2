/*
    Operator class is for :
        1. store K_parameter, 
        2. compute K, dK.
*/

#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class Operator {
protected:
    VectorXd parameter_K;
    int n_params; // how many parameters
    bool use_num_dK {false};

    SparseMatrix<double, 0, int> K, dK, d2K;
public:
    Operator(Rcpp::List ope_in)
    {
std::cout << "constructor of Operator" << std::endl;
        n_params = Rcpp::as<int> (ope_in["n_params"]);
        use_num_dK = Rcpp::as<bool> (ope_in["use_num_dK"]);
    };

    int get_n_params() const {return parameter_K.size(); }
    virtual VectorXd get_parameter() const {return parameter_K; } // kappa
    virtual void set_parameter(VectorXd parameter_K) {this->parameter_K = parameter_K;}

    // getter for K, dK, d2K
    SparseMatrix<double, 0, int>& getK()    {return K;}
    SparseMatrix<double, 0, int>& get_dK()  {return dK;}
    SparseMatrix<double, 0, int>& get_d2K() {return d2K;}

    // get K/dK using different parameter
    virtual SparseMatrix<double, 0, int> getK(VectorXd) const=0;
    virtual SparseMatrix<double, 0, int> get_dK(VectorXd) const=0;

    // param(pos) += eps;  getK(param);
    SparseMatrix<double, 0, int> getK(int pos, double eps) {
        VectorXd tmp = parameter_K;
        tmp(pos) += eps;
        return getK(tmp);
    }

    SparseMatrix<double, 0, int> get_dK(int pos, double eps) {
        VectorXd tmp = parameter_K;
        tmp(pos) += eps;
        return get_dK(tmp);
    }
    
};

// for 1. AR model 2. stationary Matern model
class matern_ope : public Operator {
private:
    // kappa = parameter(0)
    int alpha;
    VectorXd Cdiag; //, taus, kappas;
    SparseMatrix<double, 0, int> G, C;

public:
    matern_ope(Rcpp::List ope_in) 
    :   Operator    (ope_in),
        G           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]) ),
        C           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]) )
    {
        alpha = Rcpp::as<int> (ope_in["alpha"]);
        Cdiag = C.diagonal();
        
        // init parameter
            parameter_K.resize(1);
            parameter_K(0) = Rcpp::as<double> (ope_in["kappa"]);        
        
        int n_dim = G.rows();

        set_parameter(parameter_K);
        K  = getK(parameter_K);
        dK = get_dK(parameter_K);
        d2K = 0 * C;
    }
    
    void set_parameter(VectorXd parameter_K) {
        assert (parameter_K.size() == 1);
        this->parameter_K = parameter_K;

        K = getK(parameter_K);
        dK = get_dK(parameter_K);

        if (use_num_dK) {
            update_num_dK();
        }
    }

    SparseMatrix<double> getK(VectorXd params) const {
        int n_dim = G.rows();
        
        VectorXd kappas = VectorXd::Constant(n_dim, params(0));
        SparseMatrix<double> K_a (n_dim, n_dim);
        SparseMatrix<double> KCK (n_dim, n_dim);
            KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();
        
        if (alpha==2) {
            // K_a = T (G + KCK) C^(-1/2)
            // Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
            K_a = Cdiag.cwiseSqrt().cwiseInverse().asDiagonal() * (G + KCK);
        } else if (alpha==4) {      
            // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2)
            // Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
            K_a = Cdiag.cwiseSqrt().cwiseInverse().asDiagonal() * (G + KCK) * Cdiag.cwiseInverse().asDiagonal() *
            (G + KCK);
        } else {
            throw("alpha not equal to 2 or 4 is not implemented");
        }

        return K_a;
    }

    SparseMatrix<double> get_dK(VectorXd params) const {
        int n_dim = G.rows();
        SparseMatrix<double> dK (n_dim, n_dim);
        double kappa = params(0);
        if (alpha==2) 
            dK = 2*kappa*C.cwiseSqrt();
        else if (alpha==4)
            dK = 4*kappa*C.cwiseSqrt() * G + 4* pow(kappa, 3) * C.cwiseSqrt();
        else 
            throw("alpha != 2 or 4");
        return dK;
    }

    // compute numerical dK
    void update_num_dK() {
        double kappa = parameter_K(0);
        double eps = 0.01;
        SparseMatrix<double> K_add_eps = (kappa + eps) * C + G;
        dK = (K_add_eps - K) / eps;
    }
};


// Q: how to unify stationary and nonstationary case?
// for non stationary Matern model
class nonstationaryGC : public Operator {
private:
    int alpha; 
    VectorXd Cdiag; //, taus, kappas;
    SparseMatrix<double, 0, int> G;
    MatrixXd Bkappa;
public:
    nonstationaryGC(Rcpp::List ope_in) 
    :   Operator    ( ope_in)
        // alpha       ( Rcpp::as<int> (ope_in["alpha"])),
        // G           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]) )
        // Bkappa      ( Rcpp::as<MatrixXd> (ope_in["B.kappa"]) )
    {
        alpha = Rcpp::as<int> (ope_in["alpha"]);
        G = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]);
        Bkappa = Rcpp::as<MatrixXd> (ope_in["B.kappa"]);

        VectorXd theta_kappa = ope_in["init_operator"];
        SparseMatrix<double> C = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]);
        Cdiag = C.diagonal();
        set_parameter(theta_kappa);

std::cout << "Finish constructor of nonGC" << std::endl;
    }

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
            K_a = Cdiag.cwiseSqrt().cwiseInverse().asDiagonal() * (G + KCK);
        } else if (alpha==4) {      
            // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2)
            // Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
            K_a = Cdiag.cwiseSqrt().cwiseInverse().asDiagonal() * (G + KCK) * Cdiag.cwiseInverse().asDiagonal() *
            (G + KCK);
        } else {
            throw("alpha not equal to 2 or 4 is not implemented");
        }
        
        return K_a;
    }

    SparseMatrix<double> get_dK(VectorXd params) const {
        return G;        
    }

};

#endif