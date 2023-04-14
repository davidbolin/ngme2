#include "var.h"

void Var::sample_cond_V(const VectorXd& a_inc_vec, const VectorXd& b_inc_vec, int dim, bool same) {
    if (fix_V) return;
    if (noise_type == "normal") return;

    VectorXd p_vec(n), a_vec(n), b_vec(n);
    if (noise_type == "gal") {
        prevV = V;
        p_vec = (nu * h) - VectorXd::Constant(n, 0.5 * dim);
        a_vec = VectorXd::Constant(n, 2 * nu) + a_inc_vec;
        b_vec = b_inc_vec;
    } else { // nig and normal+nig
        prevV = V;
        p_vec = VectorXd::Constant(n, -0.5 - 0.5 * dim);
        a_vec = VectorXd::Constant(n, nu) + a_inc_vec;
        b_vec = VectorXd::Constant(n, nu).cwiseProduct(h).cwiseProduct(h) + b_inc_vec;
    }

    // update V
    if (same) {
        double singleV = rGIG_cpp(p_vec(0), a_vec(0), b_vec(0), var_rng());
        V = VectorXd::Constant(n, singleV);
    } else {
        V = rGIG_cpp(p_vec, a_vec, b_vec, var_rng());
    }
}

double Var::grad_log_nu() const {
    if (fix_nu || noise_type == "normal") return 0;

    // grad of log nu
    double grad = 0;
    if (noise_type == "gal") {
        for(int i=0; i < n; i++) {
            double nu_hi = nu * h[i];
            //digamma(0.1) = digamma(1.1) - 1/0.1;
            if(nu_hi > 1){
                grad -=  nu_hi * R::digamma(nu_hi);
            } else {
                grad -=  nu_hi * R::digamma(nu_hi + 1) - 1.;
            }
        grad += nu_hi * (1 - log(1/nu) + log(V(i))) - nu * V(i);
// std::cout << "grad in var = " << grad << std::endl;
        }
        grad = - grad / n;
    } else { // type == nig or normal+nig
        // df/dnu = 0.5 (2h + 1/nu - h^2/V - V)
        // df/d(log nu) = df/dnu * nu
        VectorXd tmp = 0.5 * (2*h + VectorXd::Constant(n, 1/nu)
            - h.cwiseProduct(h).cwiseQuotient(V) - V);
        double grad_nu = tmp.mean() * nu;

        VectorXd tmp2 = 0.5 * (2*h + VectorXd::Constant(n, 1/nu)
            - h.cwiseProduct(h).cwiseQuotient(prevV) - prevV);
        double grad_nu2 = tmp2.mean();

        grad = grad_nu * nu;
        double hess_nu = -0.5 * pow(nu, -2);
        // hess of log nu
        double hess = nu * grad_nu2 + nu * nu * hess_nu;

      if (hessian)
        grad = grad / hess;     // use hessian
      else
        grad = - grad / n;
// std::cout << " grad = " << grad << std::endl;
    }

    return grad;
}