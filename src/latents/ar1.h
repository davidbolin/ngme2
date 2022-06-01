#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;
class AR : public Latent {
    
public:
    AR(Rcpp::List ar1_in) 
    : Latent(ar1_in)
    {
        // read init kappa
        Rcpp::List init_value = Rcpp::as<Rcpp::List> (ar1_in["init_value"]);
        double kappa = Rcpp::as<double>  (init_value["kappa"]);

        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        ope = new GC(ope_in, kappa);
        
        // Init K and Q
        SparseMatrix<double> K = getK();
        SparseMatrix<double> Q = K.transpose() * K;
        
        solver_K.init(n_mesh, 0,0,0);
        solver_K.analyze(K);
        compute_trace();

        // Init Q
        solver_Q.init(n_mesh, 0,0,0);
        solver_Q.analyze(Q);
    }
    
    double th2k(double th) const {
        return (-1 + 2*exp(th) / (1+exp(th)));
    }
    double k2th(double k) const {
        return (log((-1-k)/(-1+k)));
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
        
        double a = ope->getKappa();
        double th = k2th(a);

        double da  = 2 * (exp(th) / pow(1+exp(th), 2));
        double d2a = 2 * (exp(th) * (-1+exp(th)) / pow(1+exp(th), 3));

        double ret = 0;
        if (numer_grad) {
            // 1. numerical gradient

            if (!use_precond) {
                double grad = (function_kappa(eps) - function_kappa(0)) / eps;
                ret = - grad * da / n_mesh;
            } else {
                double f1 = function_kappa(-eps);
                double f2 = function_kappa(0);
                double f3 = function_kappa(+eps);

                double hess = (f1 + f3 - 2*f2) / pow(eps, 2);
                double grad = (f3 - f2) / eps;
                ret = (grad * da) / (hess * da * da + grad * d2a);
            }
        } else { 
            // 2. analytical gradient and numerical hessian
            double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V) * mu);
            double grad = trace - tmp;

            if (!use_precond) {
                ret = - grad * da / n_mesh;
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

                ret = (grad * da) / (hess * da * da + grad_eps * d2a);
            }
        }

        return ret;
    }
    
    Rcpp::List get_estimates() const {
        return Rcpp::List::create(
            Rcpp::Named("alpha") = ope->getKappa(),
            Rcpp::Named("mu")    = mu,
            Rcpp::Named("sigma") = sigma,
            Rcpp::Named("var")   = var->get_var()
        );
    }
};
