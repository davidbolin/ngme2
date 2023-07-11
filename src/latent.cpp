#include "latent.h"
#include "num_diff.h"

// K is V_size * W_size matrix
Latent::Latent(const Rcpp::List& model_list, unsigned long seed) :
    latent_rng    (seed),
    model_type    (Rcpp::as<string>     (model_list["model"])),
    noise_type    (Rcpp::as<vector<string>>     (model_list["noise_type"])),
    debug         (Rcpp::as<bool>       (model_list["debug"])),
    W_size        (Rcpp::as<int>        (model_list["W_size"])),
    V_size        (Rcpp::as<int>        (model_list["V_size"])),
    n_params      (Rcpp::as<int>        (model_list["n_params"])),

    // operator
    ope           (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    ope_add_eps   (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    ope_precond   (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    h             (ope->get_h()),
    theta_K       (Rcpp::as<VectorXd>  (model_list["theta_K"])),
    n_theta_K     (theta_K.size()),
    symmetricK    (ope->is_symmetric()),
    zero_trace    (ope->is_zero_trace()),

    // K             (V_size, W_size),
    // dK            (n_theta_K),
    trace         (n_theta_K, 0.0),
    eps           (0.001),

    W             (W_size),
    prevW         (W_size),
    V             (h),
    prevV         (h),
    A             (Rcpp::as< SparseMatrix<double,0,int> >   (model_list["A"])),

    p_vec         (V_size),
    a_vec         (V_size),
    b_vec         (V_size),
    nu            (noise_type.size())
{
if (debug) std::cout << "begin constructor of latent" << std::endl;
    assert(W_size == V_size);
    // read the control variable
    Rcpp::List control_f = Rcpp::as<Rcpp::List> (model_list["control"]);
        use_precond     = Rcpp::as<bool>        (control_f["use_precond"] );
        numer_grad      = Rcpp::as<bool>        (control_f["numer_grad"]) ;
        eps             = Rcpp::as<double>      (control_f["eps"]) ;

    // construct from ngme_noise
    Rcpp::List noise_in = Rcpp::as<Rcpp::List> (model_list["noise"]);
        fix_flag[latent_fix_theta_mu]     = Rcpp::as<bool>  (noise_in["fix_theta_mu"]);
        fix_flag[latent_fix_theta_sigma]  = Rcpp::as<bool>  (noise_in["fix_theta_sigma"]);
        fix_flag[latent_fix_V]  = Rcpp::as<bool> (noise_in["fix_V"]);
        fix_flag[latent_fix_nu] = Rcpp::as<bool> (noise_in["fix_nu"]);
        single_V = Rcpp::as<bool> (noise_in["single_V"]);

        B_mu     = Rcpp::as< MatrixXd >    (noise_in["B_mu"]);
        B_sigma  = Rcpp::as< MatrixXd >    (noise_in["B_sigma"]);
        if (noise_type[0] == "normal_nig") B_sigma_normal  = Rcpp::as< MatrixXd >    (noise_in["B_sigma_normal"]);

        n_theta_mu    =   (B_mu.cols());
        n_theta_sigma =   (B_sigma.cols());
        if (noise_type[0] == "normal_nig") n_theta_sigma_normal =   (B_sigma_normal.cols());

        theta_mu = Rcpp::as< VectorXd >    (noise_in["theta_mu"]);
        theta_sigma = Rcpp::as< VectorXd > (noise_in["theta_sigma"]);
        if (noise_type[0] == "normal_nig") theta_sigma_normal = Rcpp::as< VectorXd > (noise_in["theta_sigma_normal"]);

        n_nu = Rcpp::as<int> (noise_in["n_nu"]);
        nu.resize(n_nu);
        if (n_nu > 0) nu = Rcpp::as<VectorXd> (noise_in["nu"]);

        share_V = Rcpp::as<bool> (noise_in["share_V"]);
        single_V = Rcpp::as<bool> (noise_in["single_V"]);

    // init W
    if (model_list["W"] != R_NilValue) {
        W = Rcpp::as< VectorXd > (model_list["W"]);
    } else {
        W = VectorXd::Zero(W_size);
    }
    prevW = W;
    fix_flag[latent_fix_W] = Rcpp::as<bool> (model_list["fix_W"]); // fixW

    // init K, Q
    ope->update_K(theta_K);
    if (V_size == W_size) {
        if (!symmetricK) {
            lu_solver_K.init(W_size, 0,0,0);
            lu_solver_K.analyze(getK());
        } else {
            chol_solver_K.init(W_size,0,0,0);
            chol_solver_K.analyze(getK());
        }
    }
    SparseMatrix<double> Q = getK().transpose() * getK();
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);

    // build mu, sigma, compute trace, ...
    update_each_iter(true);

    // Initialize V
    if (!Rf_isNull(noise_in["V"])) {
        V = Rcpp::as< VectorXd > (noise_in["V"]);
        prevV = V;
    } else {
        // V=h at init.
        // to-fix (for normal noise case)
        // sample_uncond_V(); // sample_uncond_V();
    }

if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
    VectorXd grad = VectorXd::Zero(n_theta_mu);
    if (n_theta_mu == 0 || fix_flag[latent_fix_theta_mu]) return grad;
    if ((n_nu == 1 && noise_type[0] == "normal") ||
        (n_nu == 2 && noise_type[0] == "normal" && noise_type[1] == "normal"))
        return grad;

    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

    // compute gradient g with V
    for (int l=0; l < n_theta_mu; l++) {
        grad(l) = (V-h).cwiseProduct(B_mu.col(l).cwiseQuotient(SV)).dot(getK()*W - mu.cwiseProduct(V-h));
    }

    return -grad;

    // // compute numerical hessian with prevV
    // VectorXd eps = VectorXd::Constant(V_size, 0.001);
    // VectorXd num_h = VectorXd::Zero(n_theta_mu);
    // for (int l=0; l < n_theta_mu; l++) {
    //     double g_o = (prevV-h).cwiseProduct(B_mu.col(l).cwiseQuotient(prevSV)).dot(getK()*W - (mu).cwiseProduct(prevV-h));
    //     double g_eps = (prevV-h).cwiseProduct(B_mu.col(l).cwiseQuotient(prevSV)).dot(getK()*W - (mu + eps).cwiseProduct(prevV-h));
    //     num_h(l) = (g_eps - g_o) / eps(0);
    // }

    // if (V_size < 10)
    //     return - grad / sqrt(W_size);
    // else
    //     return - grad / W_size;

    // double ana_hess = -(prevV-h).cwiseQuotient(prevSV).dot(prevV-h);
    // num_hess = (grad_eps - grad) / eps
}

// return the gradient wrt. theta, theta=log(sigma)
inline VectorXd Latent::grad_theta_sigma() {
    VectorXd grad = VectorXd::Zero(n_theta_sigma);
    if (fix_flag[latent_fix_theta_sigma]) return grad;

    // std::cout << "W = " << W << std::endl;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

    // (KW-mu(V-h)) / (sigma^2 V)
    VectorXd tmp = (getK()*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseQuotient(SV);

    // grad = Bi(tmp * sigma ^ -2 - 1)
    grad = B_sigma.transpose() * (tmp - VectorXd::Ones(V_size));

    return -grad;
}

inline VectorXd Latent::grad_theta_sigma_normal() {
    VectorXd V = VectorXd::Ones(V_size);
    VectorXd grad = VectorXd::Zero(n_theta_sigma_normal);

    // tmp = (KW - mu(V-h))^2 / V
    VectorXd tmp = (getK()*W).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
    // grad = Bi(tmp * sigma_normal ^ -2 - 1)
    VectorXd tmp1 = tmp.cwiseProduct(sigma_normal.array().pow(-2).matrix()) - VectorXd::Ones(V_size);
    grad += B_sigma_normal.transpose() * tmp1;

    // return - 1.0 / V_size * grad;
    return -grad;
}

VectorXd Latent::grad_theta_nu() {
    VectorXd grad = VectorXd::Zero(n_nu);

    if (n_nu == 1)
        grad(0) = NoiseUtil::grad_theta_nu(noise_type[0], nu[0], V, prevV, h, single_V);
    else if (n_nu == 2) {
        // for bivaraite case
        int n = V_size / 2;
        grad(0) = NoiseUtil::grad_theta_nu(noise_type[0], nu[0], V.segment(0, n), prevV.segment(0, n), h.segment(0, n), single_V);
        if (share_V)
            grad(1) = grad(0);
        else
            grad(1) = NoiseUtil::grad_theta_nu(noise_type[1], nu[1], V.segment(n, n), prevV.segment(n, n), h.segment(n, n), single_V);
    }
// std::cout << "g_nu = " << grad << std::endl;
    return grad;
}

VectorXd Latent::grad_theta_K() {
// std::cout << "K = " << K << std::endl;
    VectorXd grad = VectorXd::Zero(n_theta_K);
    if (fix_flag[latent_fix_theta_K]) return grad;

    if (numer_grad) {
        double val = logd_W_given_V(getK(), mu, sigma, V);
        for (int i=0; i < n_theta_K; i++) {
            VectorXd tmp = theta_K;
            tmp(i) += eps;
            ope_add_eps->update_K(tmp);
            SparseMatrix<double> K_add_eps = ope_add_eps->getK();
            double val_add_eps = logd_W_given_V(K_add_eps, mu, sigma, V);
            grad(i) = (val_add_eps - val) / eps;
        }
        // update using num_g
        // grad = num_g()
    } else {
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

        VectorXd tmp = getK() * W - mu.cwiseProduct(V-h);
        for (int j=0; j < n_theta_K; j++) {
            grad(j) = trace[j] - (ope->get_dK()[j] * W).cwiseProduct(SV.cwiseInverse()).dot(tmp);
        }
    }
// std::cout << "grad_K = " << grad.transpose() << std::endl;
    return -grad;
}

Rcpp::List Latent::output() const {
    return  Rcpp::List::create(
        Rcpp::Named("model")        = model_type,
        Rcpp::Named("noise_type")   = noise_type,
        Rcpp::Named("theta_K")      = theta_K, // same parameterization as input
        Rcpp::Named("theta_mu")     = theta_mu,
        Rcpp::Named("theta_sigma")  = theta_sigma,
        Rcpp::Named("theta_sigma_normal")  = theta_sigma_normal,
        Rcpp::Named("nu")           = nu,
        Rcpp::Named("V")            = V,
        Rcpp::Named("W")            = W
    );
}

const VectorXd Latent::get_parameter() {
if (debug) std::cout << "Start latent get parameter"<< std::endl;
    VectorXd parameter = VectorXd::Zero(n_params);

    parameter.segment(0, n_theta_K)                         = theta_K;
    parameter.segment(n_theta_K, n_theta_mu)                = theta_mu;
    parameter.segment(n_theta_K+n_theta_mu, n_theta_sigma)  = theta_sigma;
    parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_nu) = nu.array().log();
    if (noise_type[0] == "normal_nig")
        parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_nu, n_theta_sigma_normal) = theta_sigma_normal;
// if (debug) std::cout << "End latent get parameter"<< std::endl;
if (debug) std::cout << "parameter= " << parameter << std::endl;

    return parameter;
}

const VectorXd Latent::get_grad() {
// if (debug) std::cout << "Start latent get grad"<< std::endl;
    VectorXd grad = VectorXd::Zero(n_params);

    // compute gradient of each parameter
    grad.segment(0, n_theta_K) = grad_theta_K();
    grad.segment(n_theta_K, n_theta_mu) = grad_theta_mu();
    grad.segment(n_theta_K+n_theta_mu, n_theta_sigma) = grad_theta_sigma();
    grad.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_nu)  = grad_theta_nu();
    if (noise_type[0] == "normal_nig")
        grad.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_nu, n_theta_sigma_normal) = grad_theta_sigma_normal();

// if (debug) std::cout << "finish latent gradient"<< std::endl;
    return grad;
}

void Latent::set_parameter(const VectorXd& theta) {
// if (debug) std::cout << "Start latent set parameter"<< std::endl;

    // nig, gal and normal+nig
    theta_K = theta.segment(0, n_theta_K);
    theta_mu = theta.segment(n_theta_K, n_theta_mu);
    theta_sigma = theta.segment(n_theta_K+n_theta_mu, n_theta_sigma);
    nu = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma, n_nu).array().exp();
    if (noise_type[0] == "normal_nig") {
        theta_sigma_normal = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_nu, n_theta_sigma_normal);
    }
    update_each_iter();
}

void Latent::sample_cond_V() {
    if (n_nu == 0 || fix_flag[latent_fix_V]) return;
    prevV = V;

    // update b_inc (p,a_inc already built)
    b_inc = (getK() * W + mu.cwiseProduct(h)).cwiseQuotient(sigma).array().pow(2);

    int n = V_size / n_nu; // n is equal to h_size (of each mesh)
    // sample conditional V
    if (single_V) {
        if (share_V) {
            // type-G1 model (n_nu == 2, but keep same nu)
            double v1 = rGIG_cpp(-0.5-n, nu[0]+a_inc.dot(h),nu[0]+b_inc.dot(h.cwiseInverse()), latent_rng());
            V = v1 * h;
        } else {
            // type-G2 (n_nu==2) and also univariate single noise (n_nu==1)
            for (int i=0; i < n_nu; i++) {
                double v1 = rGIG_cpp(
                    -0.5 - n*1.0/n_nu,
                    nu[i] + a_inc.segment(i*n, n).dot(h.segment(i*n, n)),
                    nu[i] + b_inc.segment(i*n, n).dot(h.segment(i*n, n).cwiseInverse()),
                    latent_rng()
                );
                V.segment(i*n, n) = v1 * h.segment(i*n, n);
            }
        }
    } else {
        if (share_V) {
            // type-G3 (n_nu == 2, but keep same nu)
            // a_vec = rep(nu, n), b_vec = nu * h
            V.head(n) = rGIG_cpp(
                VectorXd::Constant(n, -1.5),
                a_vec.head(n) + a_inc.head(n) + a_inc.tail(n),
                b_vec.head(n) + b_inc.head(n) + b_inc.tail(n),
                latent_rng());
            V.tail(n) = V.head(n);
        } else {
            // type-G4 (n_nu==2) and also univariate case general noise (n_nu==1)
            for (int i=0; i < n_nu; i++) {
                if (noise_type[i] == "normal") continue;

                // sample conditional V
                V.segment(i*n, n) = rGIG_cpp((p_vec+p_inc).segment(i*n, n), (a_vec+a_inc).segment(i*n, n), (b_vec+b_inc).segment(i*n, n), latent_rng());
            }
        }

    }
}

void Latent::sample_uncond_V() {
    if (n_nu == 0 || fix_flag[latent_fix_V]) return;
    prevV = V;
    int n = V_size / n_nu;

    // same logic as in simulation.R
    if (single_V) {
        for (int i=0; i < n_nu; i++) {
            double v;
            if (noise_type[i] == "nig" || noise_type[i] == "normal_nig")
                v = rGIG_cpp(-0.5, nu[i], nu[i], latent_rng());
            else if (noise_type[i] == "gal")
                v = rGIG_cpp(nu[i], 2*nu[i], 1e-14, latent_rng());
            V.segment(i*n, n) = v * h.segment(i*n, n);
        }
    } else {
        for (int i=0; i < n_nu; i++) {
            // sample unconditional V
            V.segment(i*n, n) = rGIG_cpp(p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), latent_rng());
        }
    }

    // keep both V same for bivariate noise with share_V case
    if (n_nu == 2 && share_V) {
        V.segment(n, n) = V.segment(0, n);
    }
}

// at init, and after each set parameter
void Latent::update_each_iter(bool init) {
    if (!init) ope->update_K(theta_K);

    mu = B_mu * theta_mu;
    sigma = (B_sigma * theta_sigma).array().exp();
    if (noise_type[0] == "normal_nig")
        sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

    // update p,a,b, depend on nu, h
    if (n_nu > 0) {
        int n = V_size / n_nu;
        for (int i=0; i < n_nu; i++) {
            if (noise_type[i] == "normal") continue;
            // update p_vec, a_vec, b_vec
            NoiseUtil::update_gig(noise_type[i], nu(i),
            p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), h.segment(i*n, n), single_V);
            V.segment(i*n, n) = rGIG_cpp(p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), latent_rng());

            // update p_inc, a_inc
            p_inc = VectorXd::Constant(V_size, -0.5 * dim);
            a_inc = mu.cwiseQuotient(sigma).array().pow(2);
        }
    }

    if (!numer_grad) {
        ope->update_dK(theta_K);
        // trace[i] = tr(K^-1 dK[i])
        if (!zero_trace) {
            if (!symmetricK) {
                if (W_size > 10) {
                    lu_solver_K.computeKTK(getK());
                    for (int i=0; i < n_theta_K; i++)
                        trace[i] = lu_solver_K.trace(ope->get_dK()[i]);
                } else {
                    // for random effect case (usually small dimension)
                    for (int i=0; i < n_theta_K; i++) {
                        if (getK().toDense().isLowerTriangular() && abs(ope->get_dK()[i].diagonal().sum()) < 0.001)
                            trace[i] = 0;
                        else
                            trace[i] = getK().toDense().ldlt().solve(ope->get_dK()[i].toDense()).diagonal().sum();
                    }
                }
            } else {
                chol_solver_K.compute(getK());
                for (int i=0; i < n_theta_K; i++)
                    trace[i] = chol_solver_K.trace(ope->get_dK()[i]);
            }
        }
    }
}

// for compute hessian
double Latent::log_density(const VectorXd& parameter, bool precond_K) {
    double logd_W = 0;
    double logd_V = 0;

    if (precond_K) {
        // 1. pi(W|V)
        VectorXd theta_K = parameter.head(n_theta_K);
        VectorXd theta_mu = parameter.segment(n_theta_K, n_theta_mu);
        VectorXd theta_sigma = parameter.segment(n_theta_K + n_theta_mu, n_theta_sigma);
        VectorXd mu = B_mu * theta_mu;
        VectorXd sigma = (B_sigma * theta_sigma).array().exp();

        ope_precond->update_K(theta_K);
        SparseMatrix<double> K = ope_precond->getK();
        logd_W = logd_W_given_V(K, mu, sigma, prevV);
    } else {
        VectorXd theta_mu = parameter.head(n_theta_mu);
        VectorXd theta_sigma = parameter.segment(n_theta_mu, n_theta_sigma);
        VectorXd mu = B_mu * theta_mu;
        VectorXd sigma = (B_sigma * theta_sigma).array().exp();

        logd_W = logd_KW_given_V(mu, sigma, prevV);
    }

    // 2. pi(V)
    if ((n_nu == 1 && noise_type[0] != "normal") ||
        (n_nu == 2 && noise_type[0] != "normal" && noise_type[1] != "normal")
    ) {
        nu = parameter.tail(n_nu);
        // compute logd_V
        if (n_nu == 1) {
            logd_V = NoiseUtil::log_density(noise_type[0], V, h, nu[0], FALSE);
        } else if (n_nu == 2) {
            int n = V_size / n_nu;
            for (int i=0; i < n_nu; i++) {
                logd_V += NoiseUtil::log_density(noise_type[i], V.segment(i*n, n), h.segment(i*n, n), nu[i], FALSE);
            }
        }
    }
    return -(logd_W + logd_V);
}

// log density of W|V
double Latent::logd_W_given_V(const SparseMatrix<double>& K, const VectorXd& mu, const VectorXd& sigma, const VectorXd& V) {
    double l = 0;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

    VectorXd tmp = K * W - mu.cwiseProduct(V-h);
    if (K.rows() < 5) {
        MatrixXd Kd = K.toDense();
        l = log(Kd.diagonal().prod()) - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    } else {
        // if (!symmetricK) {
        SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
        solver_Q.compute(Q);
        l = 0.5 * solver_Q.logdet() - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        // } else {
        //     // check!!
        //     chol_solver_K.compute(K);
        //     l = chol_solver_K.logdet() - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        // }
    }
    return l;
}

// log density of KW|V
double Latent::logd_KW_given_V(const VectorXd& mu, const VectorXd& sigma, const VectorXd& V) {
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    SparseMatrix<double> K = getK();
    VectorXd tmp = K * W - mu.cwiseProduct(V-h);

    double lhs = SV.array().log().sum();
    double rhs = tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
// std::cout << "lhs: " << lhs << std::endl;
// std::cout << "rhs: " << rhs << std::endl;
    return -0.5 * (lhs + rhs);
}

// Numerical hessian
MatrixXd Latent::precond(bool precond_K) {
    VectorXd parameter (n_params);
    parameter << theta_K, theta_mu, theta_sigma, nu;

    VectorXd v (parameter);
    // update v if we don't need to compute hessian for K
    if (!precond_K) {
        v.resize(n_params - n_theta_K);
        v = parameter.tail(n_params - n_theta_K);
    }

    int n = v.size();
    MatrixXd num_hess(n, n);
	double original_val = log_density(v, precond_K);

	// compute f_v = log_density(v + precond_eps * e_i)
	VectorXd f_v (n);
	for (int i=0; i < n; i++) {
		VectorXd tmp_v = v; tmp_v(i) += precond_eps;
		f_v(i) = log_density(tmp_v, precond_K);
	}

	// compute H_ij = d2 f / dxi dxj
	for (int i=0; i < n; i++) {
		for (int j=0; j <= i; j++) {
			VectorXd tmp_vij = v; tmp_vij(i) += precond_eps; tmp_vij(j) += precond_eps;
			double f_vij = log_density(tmp_vij, precond_K);
			num_hess(i, j) = (f_vij - f_v(i) - f_v(j) + original_val) / (precond_eps * precond_eps);
		}
	}

	// fill in the lower triangular part
	for (int i=0; i < n; i++) {
		for (int j=0; j < i; j++) {
			num_hess(j, i) = num_hess(i, j);
		}
	}
// std::cout << "num_hess: " << std::endl << num_hess << std::endl;

    if (precond_K) {
        return num_hess;
    } else {
        MatrixXd precond_full = MatrixXd::Zero(n_params, n_params);
        precond_full.topLeftCorner(n_theta_K, n_theta_K) = VectorXd::Constant(n_theta_K, V_size).asDiagonal();
        precond_full.bottomRightCorner(n_params - n_theta_K, n_params - n_theta_K) = num_hess;
        return precond_full;
    }
}