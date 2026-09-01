#include "noise.h"
#include <limits>

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
    } else if (noise_type == "t" || noise_type == "skew_t") {
        // t-distribution uses inverse gamma for auxiliary variable V
        // V_i ~ InverseGamma(nu/2, nu/2), i.e., GIG(-nu/2, 0, nu)
        p = VectorXd::Constant(n, -nu(0)/2);
        a = VectorXd::Constant(n, 1e-14);
        b = VectorXd::Constant(n, nu(0));
    }
}

// compute dlog pi(V) / dtheta_nu
// Note on Rao-Blackwellising this gradient over V.
//
// It is exactly available: the gradient is affine in V and 1/V, and V | W is
// GIG(-1, a, b), so E[.|W] is closed form and cuts this gradient's variance
// about tenfold. It was implemented, measured and removed, because it does not
// pay for itself. Raising n_gibbs_samples buys more per unit of compute.
//
// Two traps if this is ever revisited: this function returns -grad (see the
// end), and a replacement must match that; and it must condition on the a, b
// actually used to draw V, not on values rebuilt from a later nu.
VectorXd NoiseUtil::grad_theta_nu(
    const string& noise_type,
    const MatrixXd& B_nu,
    const VectorXd& nu,
    const VectorXd& V,
    const VectorXd& prevV,
    const VectorXd& h,
    double nu_lower_bound,
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

            VectorXd tmp = h - V
                + h.cwiseProduct(V.array().log().matrix())
                - h.cwiseProduct(nu.cwiseInverse().array().log().matrix())
                - h.cwiseProduct(pg);

            VectorXd jac = nu.array() - nu_lower_bound;
            grad = B_nu.transpose() * tmp.cwiseProduct(jac);
        } else if (noise_type == "t" || noise_type == "skew_t") {
            // t-distribution: df/dnu = (-1 + x - x Log[2 x] + x Log[nu] - x PolyGamma[0, nu/2])/(2 x)
            // df/d(log nu) = df/dnu * nu
            VectorXd tmp(n);
            for (int j=0; j < n; j++) {
                double nu_j = nu[j];
                double x = V[j];
                tmp(j) = (-1.0 + x - x * log(2.0 * x) + x * log(nu_j) - x * R::digamma(nu_j/2.0)) / (2.0 * x);
            }
            VectorXd jac = nu.array() - nu_lower_bound;
            grad = B_nu.transpose() * tmp.cwiseProduct(jac);
        } else {
            // type == nig or normal+nig
            // df/dnu = 0.5 (2h + 1/nu - h^2/V - V)
            // df/d(log nu) = df/dnu * nu
            VectorXd tmp = 0.5 * (2*h + nu.cwiseInverse()
                - h.cwiseProduct(h).cwiseQuotient(V) - V);

            VectorXd jac = nu.array() - nu_lower_bound;
            grad = B_nu.transpose() * tmp.cwiseProduct(jac);
        }
    } else {
        // single V case
        if (noise_type == "nig") {
            // theV ~ IG(nu, nu)
            // V_i = h_i * theV
            // theV = V/h ~ GIG(-1/2, nu, nu), i.e. inverse Gaussian with mean
            // 1 and shape nu.  d/dnu log p(theV) = 1/(2 nu) - (theV-1)^2/(2 theV).
            double theV = V(0) / h(0);
            double jac = nu(0) - nu_lower_bound;
            grad(0) = (1.0 / (2.0 * nu(0))
                       - (theV - 1.0) * (theV - 1.0) / (2.0 * theV)) * jac;
        } else if (noise_type == "gal") {
            // theV ~ Gam(nu, nu)
            throw std::runtime_error("Not implemented");
        } else if (noise_type == "t") {
            // theV ~ IG(nu/2, nu/2)
            // theV = V/h ~ InverseGamma(nu/2, nu/2), so
            // d/dnu log p(theV) = 0.5*(log(nu/2) + 1 - digamma(nu/2)
            //                          - log(theV) - 1/theV).
            double theV = V(0) / h(0);
            double nu_val = nu(0);
            double jac = nu_val - nu_lower_bound;
            grad(0) = 0.5 * (std::log(nu_val / 2.0) + 1.0
                             - R::digamma(nu_val / 2.0)
                             - std::log(theV) - 1.0 / theV) * jac;
        }
    }
    return -grad;
}

// Analytic Hessian for theta_nu under NIG (and normal_nig) prior on V.
// Returns the negative Hessian, consistent with grad_theta_nu's sign.
// H = - B^T diag(nu ⊙ c) B, with c_i = h_i - h_i^2/(2 V_i) - V_i/2.
MatrixXd NoiseUtil::hess_theta_nu(
    const string& noise_type,
    const MatrixXd& B_nu,
    const VectorXd& nu,
    const VectorXd& V,
    const VectorXd& h,
    double nu_lower_bound,
    bool single_V
) {
    const int n_theta = B_nu.cols();
    MatrixXd H = MatrixXd::Zero(n_theta, n_theta);

    // Guard conditions consistent with grad implementation
    if (nu.size() == 0 || B_nu.cols() == 0) return H;
    if (nu.mean() > 10000) return H;
    if (noise_type == "normal") return H; // no contribution

    if (!single_V) {
        if (noise_type == "nig" || noise_type == "normal_nig") {
            // c = h - (h^2)/(2V) - V/2
            VectorXd c = h - 0.5 * h.cwiseProduct(h).cwiseQuotient(V) - 0.5 * V;
            VectorXd d = (nu.array() - nu_lower_bound).matrix().cwiseProduct(c); // diag entries
            H = - (B_nu.transpose() * d.asDiagonal() * B_nu);
        } else if (noise_type == "gal") {
            // GAL: V ~ Gamma(h*nu, scale=1/nu)
            // g_i(nu_i) = h_i log v_i - v_i + h_i (log nu_i + 1) - h_i psi(h_i nu_i)
            // d_i = nu_i g_i(nu_i) + h_i nu_i - h_i^2 nu_i^2 psi^{(1)}(h_i nu_i)
            const int n = V.size();
            VectorXd logV = V.array().log().matrix();
            VectorXd dig(n), trig(n);
            for (int i = 0; i < n; ++i) {
                double hx = h(i) * nu(i);
                dig(i)  = R::digamma(hx);
                // trigamma = derivative of digamma = psigamma(x, 1)
                trig(i) = R::psigamma(hx, 1.0);
            }
            VectorXd g = h.cwiseProduct(logV) - V
                       + h.cwiseProduct(nu.array().log().matrix() + VectorXd::Ones(n))
                       - h.cwiseProduct(dig);
            VectorXd d = nu.cwiseProduct(g)
                       + h.cwiseProduct(nu)
                       - h.cwiseProduct(h).cwiseProduct(nu.cwiseProduct(nu)).cwiseProduct(trig);
            H = - (B_nu.transpose() * d.asDiagonal() * B_nu);
        } else {
            // GAL / t-skewt not implemented analytically here
            // Keep zero so numerical path (if used elsewhere) dominates
        }
    } else {
        // Single-V special model uses a different reparam; leave to numeric Hessian
        // to avoid mismatch with grad's special-case formula.
    }
    return H;
}

double NoiseUtil::log_density(
    const string& noise_type,
    const VectorXd& V,
    const VectorXd& h,
    const MatrixXd& B_nu,
    const VectorXd& theta_nu,
    bool single_V
) {
    if (noise_type != "nig" && noise_type != "gal" && noise_type != "t") return 0;

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
            logd += pow(beta, alpha) / R::gammafn(alpha) * pow(x, alpha-1) * exp(-beta*x);
        } else if (noise_type == "t") {
            // V_i ~ IG(nu/2, nu/2) for t-distribution auxiliary variable
            double alpha = nu(i) / 2;
            double beta = nu(i) / 2;
            logd += alpha * log(beta) - R::lgammafn(alpha) - (alpha + 1) * log(x) - beta / x;
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
