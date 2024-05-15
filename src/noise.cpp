#include "noise.h"

void NoiseUtil::update_gig(
    const string& noise_type,
    const VectorXd& nu,
    Eigen::Ref<Eigen::VectorXd> p,
    Eigen::Ref<Eigen::VectorXd> a,
    Eigen::Ref<Eigen::VectorXd> b,
    const VectorXd& h,
    bool single_V
) {
    int n = h.size();
    if (noise_type == "gal") {
        p = h.cwiseProduct(nu);
        a = 2 * nu;
        b = VectorXd::Constant(n, 1e-14);
    } else if (noise_type == "nig" || noise_type == "normal_nig") {
        p = VectorXd::Constant(n, -0.5);
        if (!single_V) {
            a = nu;
            b = a.cwiseProduct(h.cwiseProduct(h));
        } else {
            // V_tilde ~ IG(nu/h, nu h)
            a = nu.cwiseProduct(h.cwiseInverse());
            b = nu.cwiseProduct(h);
        }
    }
}

// compute dlog pi(V) / dtheta_nu
VectorXd NoiseUtil::grad_theta_nu(
    const string& noise_type,
    const MatrixXd& B_nu,
    const VectorXd& nu,
    const VectorXd& V,
    const VectorXd& prevV,
    const VectorXd& h,
    bool single_V
) {
    int n_theta_nu = B_nu.cols();
    VectorXd grad = VectorXd::Zero(n_theta_nu);

    if (nu.mean() > 10000) return grad;
    if (noise_type == "normal") return grad;

    int n = V.size();
    if (!single_V) {
        if (noise_type == "gal") {
            VectorXd pg (n);
            for (int j=0; j < n; j++) pg(j) = R::digamma(nu[j]*h[j]);
            VectorXd tmp = h - V + h.cwiseProduct(V.array().log().matrix())
                - h.cwiseProduct(nu.cwiseInverse().array().log().matrix())
                - h.cwiseProduct(pg);
            grad = B_nu.transpose() * tmp.cwiseProduct(nu);
            // for (int j=0; j < n; j++) {
            //     double nu_hj = nu[j] * h[j];
            //     double v = V(j);
            //     grad(i) -= (-v + h[j] * (1 + log(v) - log(1/nu[j]) + R::digamma(nu_hj))) * B_nu(j, i);
            //     grad(i) /= n;
            //     std::cout << "grad: " << grad(i) << std::endl;
            // }
        } else {
            // type == nig or normal+nig
            // df/dnu = 0.5 (2h + 1/nu - h^2/V - V)
            // df/d(log nu) = df/dnu * nu
            VectorXd tmp = 0.5 * (2*h + nu.cwiseInverse()
                - h.cwiseProduct(h).cwiseQuotient(V) - V);
            grad = B_nu.transpose() * tmp.cwiseProduct(nu);
        // skip the hessian
            // VectorXd tmp2 = 0.5 * (2*h + VectorXd::Constant(n, 1/nu)
            //     - h.cwiseProduct(h).cwiseQuotient(prevV) - prevV);
            // double grad_nu2 = tmp2.mean();

            // grad = grad_nu * nu;
            // double hess_nu = -0.5 * pow(nu, -2);
            // // hess of log nu
            // double hess = nu * grad_nu2 + nu * nu * hess_nu;

            // grad *= n;
        }
    } else {
        // single V case
        if (noise_type == "nig") {
            // theV ~ IG(nu, nu)
            // V_i = h_i * theV
            double theV = V(0) / h(0);
            grad(0) = - 0.1 * (nu(0) - 3*theV - nu(0)*theV*theV)/(2*theV*theV); // remove negative sign
        } else if (noise_type == "gal") {
            // theV ~ Gam(nu, nu)
            throw std::runtime_error("Not implemented");
        }
    }
    // grad /= n;
    return -grad;
}

double NoiseUtil::log_density(
    const string& noise_type,
    const VectorXd& V,
    const VectorXd& h,
    const MatrixXd& B_nu,
    const VectorXd& theta_nu,
    bool single_V
) {
    if (noise_type != "nig" && noise_type != "gal") return 0;

    VectorXd nu = (B_nu * theta_nu).array().exp();
    assert(V.size() == h.size());
    double logd=0;
    for (int i = 0; i < V.size(); i++) {
        double x = V(i);
        if (noise_type == "nig") {
            // V_i ~ IG(nu, nu h_i)
            double mu = nu(i);
            double lambda = nu(i) * h(i);
            logd += 0.5 * log(lambda / (2*Pi*pow(x, 3)))
                - lambda * (x - mu) / (2*mu*mu*x);
        } else if (noise_type == "gal") {
            // V_i ~ Gamma(alpha=h_i nu, beta=nu), see wikipedia
            double alpha = h(i) * nu(i);
            double beta = nu(i);
            logd += pow(beta, beta) / R::gammafn(alpha) * pow(x, alpha-1) * exp(-beta*x);
        }
    }

    return logd;
}

// seems off-diagonal is 0
MatrixXd NoiseUtil::precond(
    double eps,
    const string& noise_type,
    const VectorXd& V,
    const VectorXd& h,
    const MatrixXd& B_nu,
    const VectorXd& theta_nu,
    bool single_V
) {
    // (f(x+h) - 2f(x) +f(x-h)) / h^2
    // double f1 = log_density(noise_type, V, h, nu+eps, single_V);
    // double f2 = log_density(noise_type, V, h, nu, single_V);
    // double f3 = log_density(noise_type, V, h, nu-eps, single_V);
    // return - (f1 - 2*f2 + f3) / (eps*eps);

    int n = theta_nu.size();
    int n_theta_nu = B_nu.cols();

    VectorXd nu = (B_nu * theta_nu).array().exp();
    MatrixXd hessian = MatrixXd::Zero(n_theta_nu, n_theta_nu);
    double original_val = log_density(noise_type, V, h, B_nu, theta_nu, single_V);

	// compute f_v = log_density(v + eps * e_i)
	VectorXd f_v (n_theta_nu);
	for (int i=0; i < n; i++) {
		VectorXd tmp_v = theta_nu; tmp_v(i) += eps;
		f_v(i) = log_density(noise_type, V, h, B_nu, tmp_v, single_V);
	}

	// compute H_ij = d2 f / dxi dxj
	for (int i=0; i < n; i++) {
		for (int j=0; j <= i; j++) {
			VectorXd tmp_vij = theta_nu; tmp_vij(i) += eps; tmp_vij(j) += eps;
			double f_vij = log_density(noise_type, V, h, B_nu, tmp_vij, single_V);
			hessian(i, j) = (f_vij - f_v(i) - f_v(j) + original_val) / (eps * eps);
		}
	}

	// fill in the lower triangular part
	for (int i=0; i < n; i++) {
		for (int j=0; j < i; j++) {
			hessian(j, i) = hessian(i, j);
		}
	}

    return -hessian;
}