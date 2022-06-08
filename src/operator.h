#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

using Eigen::SparseMatrix;
using Eigen::VectorXd;

class Operator {
protected:
    VectorXd parameter;
    int n_params;
    bool numerical_dK {false};

    SparseMatrix<double, 0, int> K, dK, d2K;
public:
    Operator(Rcpp::List ope_in) {
        n_params = Rcpp::as<bool> (ope_in["n_params"]);
std::cout << "n_params in latent = " << n_params << std::endl;
        numerical_dK = Rcpp::as<bool> (ope_in["numerical_dK"]);
    };

    virtual VectorXd get_parameter() const {return parameter; }
    int get_n_params() const {return parameter.size(); }
    
    virtual void set_parameter(VectorXd parameter)=0;

    // getter for K, dK, d2K
    SparseMatrix<double, 0, int>& getK()    {return K;}
    SparseMatrix<double, 0, int>& get_dK()  {return dK;}
    SparseMatrix<double, 0, int>& get_d2K() {return d2K;}

    // K(kappa + eps), dK(kappa + eps)
    virtual SparseMatrix<double, 0, int> getK(double eps)   const=0;
    virtual SparseMatrix<double, 0, int> get_dK(double eps) const=0;
};

// fit for stationary AR and Matern 
class stationaryGC : public Operator {
private:
    // kappa = parameter(0)
    SparseMatrix<double, 0, int> G, C;

public:
    stationaryGC(Rcpp::List ope_in, double kappa) 
    : Operator(ope_in) 
    {
        // init parameter
            parameter.resize(1);
            parameter(0) = kappa;

        G = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]);
        C = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]);
        
        K =  kappa * C + G;
        dK = C;
        d2K = 0 * C;
    }

    void set_parameter(VectorXd parameter) {
        assert (parameter.size() == 1);
        this->parameter = parameter;

        K = parameter(0) * C + G;
        if (numerical_dK) {
            update_num_dK();
        }
    }

    SparseMatrix<double> getK(double eps) const {
        double kappa = parameter(0);
        SparseMatrix<double> K =  (kappa+eps) * C + G;
        return K;
    }

    SparseMatrix<double> get_dK(double eps) const {
        return C;
    }

    void update_num_dK() {
        double kappa = parameter(0);
        double eps = 0.01;
        SparseMatrix<double> Keps = (kappa + eps) * C + G;
        dK = (Keps - K) / eps;
    }

};

// class nonstationaryGC : public Operator {
// private:
//     // kappa = parameter(0)
//     SparseMatrix<double, 0, int> G, C;

// public:
//     stationaryGC(Rcpp::List ope_in, double kappa) 
//     : Operator(ope_in) 
//     {
//         // init parameter
//             parameter.resize(1);
//             parameter(0) = kappa;

//         G = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]);
//         C = Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]);
        
//         K =  kappa * C + G;
//         dK = C;
//         d2K = 0 * C;
//     }

//     void set_parameter(VectorXd parameter) {
//         assert (parameter.size() == 1);
//         this->parameter = parameter;

//         K = parameter(0) * C + G;
//         if (numerical_dK) {
//             update_num_dK();
//         }
//     }

//     SparseMatrix<double> getK(double eps) const {
//         double kappa = parameter(0);
//         SparseMatrix<double> K =  (kappa+eps) * C + G;
//         return K;
//     }

//     SparseMatrix<double> get_dK(double eps) const {
//         return C;
//     }

//     void update_num_dK() {
//         double kappa = parameter(0);
//         double eps = 0.01;
//         SparseMatrix<double> Keps = (kappa + eps) * C + G;
//         dK = (Keps - K) / eps;
//     }

// };

#endif