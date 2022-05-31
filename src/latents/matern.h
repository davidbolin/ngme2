#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;
class Matern : public Latent {
    
public:
    Matern(Rcpp::List matern_in) 
    : Latent(matern_in)
    {
        // read init kappa
        Rcpp::List init_value = Rcpp::as<Rcpp::List> (matern_in["init_value"]);
        double kappa = Rcpp::as<double>  (init_value["kappa"]);

        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (matern_in["operator_in"]); // containing C and G
        ope = new GC(ope_in, kappa);
        
        // Init K and Q
        SparseMatrix<double> K = getK();
        SparseMatrix<double> Q = K.transpose() * K;
        
        solver_K.init(n_reg, 0,0,0);
        solver_K.analyze(K);
        compute_trace();

        // Init Q
        solver_Q.init(n_reg, 0,0,0);
        solver_Q.analyze(Q);
    }
    
    double th2k(double th) const {
        return exp(th);
    }
    double k2th(double k) const {
        return log(k);
    }

    double get_theta_kappa() const {
        double k = ope->getKappa();
        return k2th(k);
    }

    void set_theta_kappa(double v) {
        double k = th2k(v);
        setKappa(k);
    }

    // grad_kappa * dkappa/dtheta
    double grad_theta_kappa() {
        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
        VectorXd V = getV();
        VectorXd SV = pow(sigma, 2) * V;
        
        double k = ope->getKappa();
        double th = k2th(k);

        double dk  = k;
        double d2k = k;

        double ret = 0;
        if (numer_grad) {
            // 1. numerical gradient
            double kappa = ope->getKappa();

            if (!use_precond) {
                double grad = (function_kappa(eps) - function_kappa(0)) / eps;
                ret = - grad * dk / n_reg;
            } else {
                double f1 = function_kappa(-eps);
                double f2 = function_kappa(0);
                double f3 = function_kappa(+eps);

                double hess = (f1 + f3 - 2*f2) / pow(eps, 2);
                double grad = (f3 - f2) / eps;
                ret = (grad * dk) / (hess * dk * dk + grad * d2k);
            }
        } else { 
            // 2. analytical gradient and numerical hessian
            double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V) * mu);
            double grad = trace - tmp;

            if (!use_precond) {
                ret = - grad * dk / n_reg;
            } else {
                VectorXd prevV = getPrevV();
                // compute numerical hessian
                SparseMatrix<double> K2 = ope->getK(eps);
                SparseMatrix<double> dK2 = ope->get_dK(eps);

                // grad(x+eps) - grad(x) / eps
                VectorXd prevSV = pow(sigma, 2) * prevV;
                double grad2_eps = trace_eps - (dK2*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K2 * prevW + (h - prevV) * mu);
                double grad_eps  = trace - (dK*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K * prevW + (h - prevV) * mu);

                double hess = (grad2_eps - grad_eps) / eps;

                ret = (grad * dk) / (hess * dk * dk + grad_eps * d2k);
            }
        }

        return ret;
    }
    
    Rcpp::List get_estimates() const {
        return Rcpp::List::create(
            Rcpp::Named("kappa") = ope->getKappa(),
            Rcpp::Named("mu")    = mu,
            Rcpp::Named("sigma") = sigma,
            Rcpp::Named("var")   = var->get_var()
        );
    }
};
