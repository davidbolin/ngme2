// AR1 and its operator

#ifndef NGME_AR1_H
#define NGME_AR1_H

#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include "../operator.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;

// For AR1 model: parameter_K(0) = alpha (everything wrt. alpha)
class ar_operator : public Operator {
private:
    SparseMatrix<double, 0, int> G, C;

public:
    ar_operator(Rcpp::List ope_in) 
    :   Operator    (ope_in),
        G           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["G"]) ),
        C           ( Rcpp::as< SparseMatrix<double,0,int> > (ope_in["C"]) )
    {
        // init parameter
            parameter_K.resize(1);
            parameter_K(0) = Rcpp::as<double> (ope_in["alpha"]);

        K  = getK(parameter_K);
        dK = get_dK(parameter_K);
        d2K = 0 * C;
    }
        
    void set_parameter(VectorXd alpha) {
        assert (alpha.size() == 1);
        this->parameter_K = alpha;

        K = getK(alpha);
        if (use_num_dK) {
            update_num_dK();
        }
    }

    // export 
    SparseMatrix<double> getK(VectorXd alpha) const {
        assert (alpha.size() == 1);
        SparseMatrix<double> K = alpha(0) * C + G;
        return K;
    }

    SparseMatrix<double> get_dK(VectorXd alpha) const {
        return C;
    }

    // compute numerical dK
    void update_num_dK() {
        double alpha = parameter_K(0);
        double eps = 0.01;
        SparseMatrix<double> K_add_eps = (alpha + eps) * C + G;
        dK = (K_add_eps - K) / eps;
    }
};


// get_K_params, grad_K_params, set_K_params, output
class AR : public Latent {
public:
    AR(Rcpp::List latent_in) 
    : Latent(latent_in)
    {
        Rcpp::List operator_in = Rcpp::as<Rcpp::List> (latent_in["operator_in"]); // containing C and G
        
        // Init operator for ar1
        ope = new ar_operator(operator_in);
        
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
    
    // override get_the_K with change of variable
    VectorXd get_theta_K() const {
        VectorXd alpha = ope->get_parameter();
        assert (alpha.size() == 1);
        // change of variable
        double th = a2th(alpha(0));
        return VectorXd::Constant(1, th);
    } 
    // return length 1 vectorxd : grad_kappa * dkappa/dtheta 
    VectorXd grad_theta_K() {
std::cout << "begin ar1 grad_theta_kappa" << std::endl;

        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
        VectorXd V = getV();
        VectorXd SV = getSV();
        
        VectorXd params = ope->get_parameter();
            double a = params(0);
        double th = a2th(a);

        double da  = 2 * (exp(th) / pow(1+exp(th), 2));
        double d2a = 2 * (exp(th) * (-1+exp(th)) / pow(1+exp(th), 3));

        double ret = 0;
        if (numer_grad) {
            // 1. numerical gradient
            if (!use_precond) {
                // double grad = (function_K(eps) - function_K(0)) / eps;
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
            double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
            double grad = trace - tmp;

            if (!use_precond) {
                ret = - grad * da / n_mesh;
            } else {
                VectorXd prevV = getPrevV();
                // compute numerical hessian
                SparseMatrix<double> K2 = ope->getK(0, eps);
                SparseMatrix<double> dK2 = ope->get_dK(0, eps);

                // grad(x+eps) - grad(x) / eps
                VectorXd prevSV = getPrevSV();
                double grad2_eps = trace_eps - (dK2*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K2 * prevW + (h - prevV).cwiseProduct(mu));
                double grad_eps  = trace - (dK*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K * prevW + (h - prevV).cwiseProduct(mu));

                double hess = (grad2_eps - grad_eps) / eps;

                ret = (grad * da) / (hess * da * da + grad_eps * d2a);
            }
        }
        
        return VectorXd::Constant(1, ret);
    }

    void set_theta_K(VectorXd theta) {
        // change of variable
        double alpha = th2a(theta(0));
        ope->set_parameter(VectorXd::Constant(1, alpha));

        if (!numer_grad) compute_trace(); 
    }

    double th2a(double th) const {
        return (-1 + 2*exp(th) / (1+exp(th)));
    }
    
    double a2th(double k) const {
        return (log((-1-k)/(-1+k)));
    }
    
    // generating output
    Rcpp::List get_estimates() const {
        return Rcpp::List::create(
            Rcpp::Named("alpha")        = ope->get_parameter()(0),
            Rcpp::Named("theta.mu")     = theta_mu,
            Rcpp::Named("theta.sigma")  = theta_sigma,
            Rcpp::Named("theta.noise")  = var->get_var()
        );
    }
};

#endif