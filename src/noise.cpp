#include "noise.h"

void NoiseUtil::update_gig(
    const string& noise_type,
    double nu,
    Eigen::Ref<Eigen::VectorXd> p,
    Eigen::Ref<Eigen::VectorXd> a,
    Eigen::Ref<Eigen::VectorXd> b,
    const VectorXd& h,
    bool single_V
) {
    int n = h.size();
    if (noise_type == "gal") {
        p = h * nu;
        a = VectorXd::Constant(n, nu * 2);
        b = VectorXd::Constant(n, 1e-14);
    } else if (noise_type == "nig" || noise_type == "normal_nig") {
        p = VectorXd::Constant(n, -0.5);
        if (!single_V) {
            a = VectorXd::Constant(n, nu);
            b = a.cwiseProduct(h.cwiseProduct(h));
        } else {
            // V_tilde ~ IG(nu/h, nu h)
            a = nu * h.cwiseInverse();
            b = nu * h;
        }
    }
}

double NoiseUtil::grad_theta_nu(
    const string& noise_type,
    double nu,
    const VectorXd& V,
    const VectorXd& prevV,
    const VectorXd& h,
    bool single_V
) {
    if (noise_type == "normal") return 0;
    int n = V.size();
    bool use_hessian = true;
    double grad = 0;

    if (!single_V) {
        if (noise_type == "gal") {
            for (int i=0; i < n; i++) {
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
        } else {
            // type == nig or normal+nig
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

        if (use_hessian)
            grad = grad / hess;
        else
            grad = - grad / n;
        }
    } else {
        // single V case
        if (noise_type == "nig") {
            // theV ~ IG(nu, nu)
            // V_i = h_i * theV
            double theV = V(0) / h(0);
            double theV2 = V(1) / h(1);
// std::cout << "theV = " << theV << std::endl; std::cout << "theV2 = " << theV2 << std::endl;
            grad = (nu - 3*theV - nu*theV*theV)/(2*theV*theV); // remove negative sign
        } else if (noise_type == "gal") {
            // to-do
        }
    }

    return grad;
}

