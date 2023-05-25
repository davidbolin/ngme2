#include "noise.h"

void NoiseUtil::update_gig(
    const string& noise_type,
    double nu,
    Eigen::Ref<Eigen::VectorXd> p,
    Eigen::Ref<Eigen::VectorXd> a,
    Eigen::Ref<Eigen::VectorXd> b,
    const VectorXd& h
) {
    int n = h.size();
    if (noise_type == "gal") {
        p = h * nu;
        a = VectorXd::Constant(n, nu * 2);
        b = VectorXd::Constant(n, 1e-14);
    } else if (noise_type == "nig" || noise_type == "normal_nig") {
        p = VectorXd::Constant(n, -0.5);
        a = VectorXd::Constant(n, nu);
        b = a.cwiseProduct(h.cwiseProduct(h));
    }
}

double NoiseUtil::grad_theta_nu(
    const string& noise_type,
    double nu,
    const VectorXd& V,
    const VectorXd& prevV,
    const VectorXd& h
) {
    if (noise_type == "normal") return 0;
    int n = V.size();
    bool hessian = true;

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
    }
    return grad;
}

