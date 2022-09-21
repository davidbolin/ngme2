// Notice here:
// for R interface, theta_K is directly the parameter_K of the Opeartor object
// for C interface, theta_K is f(parameter_K) such that it is unbounded

// AR1 and its operator

#ifndef NGME_AR1_H
#define NGME_AR1_H

#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;

/*
    AR model:
        parameter_K(0) = alpha
        K = C * alpha + G
*/

// W_size = V_size
// get_K_params, grad_K_params, set_K_params, output
class AR : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
public:
    AR(Rcpp::List& model_list, unsigned long seed)
    : Latent(model_list, seed),
      G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
      C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"]))
    {
if (debug) std::cout << "Begin Constructor of AR1" << std::endl;
        // Init operator for ar1
        // ope = new ar_operator(model_list);

        // VectorXd parameter_K = Rcpp::as<VectorXd> (model_list["theta_K"]);
        // set_parameter(parameter_K);

        // Init K and Q
        // SparseMatrix<double> K = getK();
        K = getK(parameter_K);
        SparseMatrix<double> Q = K.transpose() * K;

        // watch out!
        if (W_size == V_size) {
            lu_solver_K.init(W_size, 0,0,0);
            lu_solver_K.analyze(K);
            compute_trace();
        }

        // Init Q
        solver_Q.init(W_size, 0,0,0);
        solver_Q.analyze(Q);
if (debug) std::cout << "End Constructor of AR1" << std::endl;
    }

    // For optimizer
    VectorXd get_unbounded_theta_K() const {
        VectorXd alpha = parameter_K;
        assert (alpha.size() == 1);
        // change of variable
        double th = a2th(alpha(0));
        return VectorXd::Constant(1, th);
    }

    // wrt. parameter_K (bounded parameter)
    SparseMatrix<double> getK(VectorXd alpha) const {
        assert (alpha.size() == 1);
        SparseMatrix<double> K = alpha(0) * C + G;
        return K;
    }

    SparseMatrix<double> get_dK(int index, VectorXd alpha) const {
        assert(index==0);
        return C;
    }

    // compute numerical dK
    void update_num_dK() {
        double alpha = parameter_K(0);
        double eps = 0.01;
        SparseMatrix<double> K_add_eps = (alpha + eps) * C + G;
        dK = (K_add_eps - K) / eps;
    }

    // return length 1 vectorxd : grad_kappa * dkappa/dtheta
    VectorXd grad_theta_K() {
        SparseMatrix<double> dK = get_dK_by_index(0);
        VectorXd V = getV();
        VectorXd SV = getSV();

        double a = parameter_K(0);
        double th = a2th(a);

        double da  = 2 * (exp(th) / pow(1+exp(th), 2));
        double d2a = 2 * (exp(th) * (-1+exp(th)) / pow(1+exp(th), 3));

        double ret = 0;
        if (numer_grad) {
            // 1. numerical gradient
            ret = numerical_grad()(0);
        } else {
            // 2. analytical gradient and numerical hessian
            double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
            double grad = trace - tmp;

            if (!use_precond) {
                ret = - grad * da / W_size;
            } else {
                VectorXd prevV = getPrevV();
                // compute numerical hessian
                SparseMatrix<double> K2 = getK_by_eps(0, eps);
                SparseMatrix<double> dK2 = get_dK_by_eps(0, 0, eps);

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

    void set_unbound_theta_K(VectorXd theta) {
        // change of variable
        double alpha = th2a(theta(0));
        parameter_K = VectorXd::Constant(1, alpha);
        K = getK(parameter_K);
        dK = get_dK(0, parameter_K);
        d2K = 0 * C;

        if (use_num_dK) {
            update_num_dK();
        }

        if (!numer_grad && (W_size == V_size)) compute_trace();
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
            Rcpp::Named("alpha")        = parameter_K(0),
            Rcpp::Named("theta.mu")     = theta_mu,
            Rcpp::Named("theta.sigma")  = theta_sigma,
            Rcpp::Named("theta.noise")  = var->get_var()
        );
    }
};

#endif