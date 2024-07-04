#include "latent.h"
#include "prior.h"
#include "num_diff.h"

// K is V_size * W_size matrix
Latent::Latent(const Rcpp::List& model_list, unsigned long seed) :
    latent_rng    (seed),
    model_type    (Rcpp::as<string>     (model_list["model"])),
    noise_type    (Rcpp::as<vector<string>>     (model_list["noise_type"])),
    debug         (Rcpp::as<bool>       (model_list["debug"])),
    n_noise       (noise_type.size()),
    W_size        (Rcpp::as<int>        (model_list["W_size"])),
    V_size        (Rcpp::as<int>        (model_list["V_size"])),
    n_params      (Rcpp::as<int>        (model_list["n_params"])),

    // operator
    ope           (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    ope_precond   (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    ope_addeps    (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    h             (ope->get_h()),
    theta_K       (Rcpp::as<VectorXd>  (model_list["theta_K"])),
    n_theta_K     (theta_K.size()),
    symmetricK    (ope->is_symmetric()),
    zero_trace    (ope->is_zero_trace()),

    // K             (V_size, W_size),
    // dK            (n_theta_K),
    trace         (n_theta_K, 0.0),
    rb_trace_K    (n_theta_K),
    eps           (0.001),
    num_dK        (n_theta_K, SparseMatrix<double>(V_size, W_size)),

    W             (W_size),
    prevW         (W_size),
    cond_W        (W_size),
    V             (h),
    prevV         (h),
    A             (Rcpp::as<SparseMatrix<double,0,int>>   (model_list["A"])),

    p_vec         (V_size),
    a_vec         (V_size),
    b_vec         (V_size)
{
if (debug) std::cout << "begin constructor of latent" << std::endl;
    assert(W_size == V_size);
    // read the control variable
    Rcpp::List control_f = Rcpp::as<Rcpp::List> (model_list["control"]);
        numer_grad      = Rcpp::as<bool>        (control_f["numer_grad"]) ;
if (control_f.containsElementNamed("improve_hessian"))
        improve_hessian = Rcpp::as<bool>        (control_f["improve_hessian"]) ;
        eps             = Rcpp::as<double>      (control_f["eps"]) ;

    // construct from ngme_noise
    fix_flag[latent_fix_theta_K]     = Rcpp::as<bool>  (model_list["fix_theta_K"]);
    Rcpp::List noise_in = Rcpp::as<Rcpp::List> (model_list["noise"]);
        fix_flag[latent_fix_theta_mu]     = Rcpp::as<bool>  (noise_in["fix_theta_mu"]);
        fix_flag[latent_fix_theta_sigma]  = Rcpp::as<bool>  (noise_in["fix_theta_sigma"]);
if (noise_in.containsElementNamed("fix_theta_nu"))
        fix_flag[latent_fix_theta_nu]     = Rcpp::as<bool> (noise_in["fix_theta_nu"]);
if (noise_in.containsElementNamed("latent_fix_theta_sigma_normal"))
        fix_flag[latent_fix_theta_sigma_normal]  = Rcpp::as<bool> (noise_in["fix_theta_sigma_normal"]);
        fix_flag[latent_fix_V]  = Rcpp::as<bool> (noise_in["fix_V"]);

        single_V = Rcpp::as<bool> (noise_in["single_V"]);

        B_mu     = Rcpp::as< MatrixXd >    (noise_in["B_mu"]);
        B_sigma  = Rcpp::as< MatrixXd >    (noise_in["B_sigma"]);
        B_nu     = Rcpp::as< MatrixXd >    (noise_in["B_nu"]);

        if (noise_type[0] == "normal_nig") B_sigma_normal  = Rcpp::as< MatrixXd >    (noise_in["B_sigma_normal"]);

        n_theta_mu    = Rcpp::as< int > (noise_in["n_theta_mu"]);
        n_theta_sigma = Rcpp::as< int > (noise_in["n_theta_sigma"]);
        n_theta_nu    = Rcpp::as< int > (noise_in["n_theta_nu"]);

        if (noise_in.containsElementNamed("nu_lower_bound"))
            nu_lower_bound  = Rcpp::as< double > (noise_in["nu_lower_bound"]);

        rb_trace_sigma = VectorXd::Zero(n_theta_sigma);
        if (noise_type[0] == "normal_nig") n_theta_sigma_normal =   (B_sigma_normal.cols());

        theta_mu    = Rcpp::as< VectorXd > (noise_in["theta_mu"]);
        theta_sigma = Rcpp::as< VectorXd > (noise_in["theta_sigma"]);
        theta_nu    = Rcpp::as< VectorXd > (noise_in["theta_nu"]);
        if (noise_type[0] == "normal_nig") theta_sigma_normal = Rcpp::as< VectorXd > (noise_in["theta_sigma_normal"]);

        share_V = Rcpp::as<bool> (noise_in["share_V"]);
        single_V = Rcpp::as<bool> (noise_in["single_V"]);

        // init priors for parameter of K and noise
        Rcpp::List prior_list = Rcpp::as<Rcpp::List> (model_list["prior_theta_K"]);
            prior_K_type  = Rcpp::as<string> (prior_list["type"]);
            prior_K_param = Rcpp::as<VectorXd> (prior_list["param"]);
        prior_list = Rcpp::as<Rcpp::List> (noise_in["prior_mu"]);
            prior_mu_type  = Rcpp::as<string> (prior_list["type"]);
            prior_mu_param = Rcpp::as<VectorXd> (prior_list["param"]);
        prior_list = Rcpp::as<Rcpp::List> (noise_in["prior_sigma"]);
            prior_sigma_type  = Rcpp::as<string> (prior_list["type"]);
            prior_sigma_param = Rcpp::as<VectorXd> (prior_list["param"]);
        prior_list = Rcpp::as<Rcpp::List> (noise_in["prior_nu"]);
            prior_nu_type  = Rcpp::as<string> (prior_list["type"]);
            prior_nu_param = Rcpp::as<VectorXd> (prior_list["param"]);

// std::cout << " here 2" << std::endl;

    // init W
    if (model_list["W"] != R_NilValue) {
        W = Rcpp::as< VectorXd > (model_list["W"]);
    } else {
        W = VectorXd::Zero(W_size);
    }
    prevW = W;
    fix_flag[latent_fix_W] = Rcpp::as<bool> (model_list["fix_W"]); // fixW

    // init K, Q, dK
    int n_trace_iter = Rcpp::as<int> (model_list["n_trace_iter"]);
    ope->update_K(theta_K);
    if (V_size == W_size) {
        if (!symmetricK) {
            // lu_solver_K.set_N(n_trace_iter);
            lu_solver_K.init(W_size, n_trace_iter, 0, 0);
            lu_solver_K.analyze(getK());
        } else {
            // chol_solver_K.set_N(n_trace_iter);
            chol_solver_K.init(W_size, n_trace_iter, 0, 0);
            chol_solver_K.analyze(getK());
        }
    }
    SparseMatrix<double> Q = getK().transpose() * getK();
    solver_Q.init(W_size, n_trace_iter,0,0);
    solver_Q.analyze(Q);

// std::cout << " here 3" << std::endl;
    // build mu, sigma, compute trace, ...
    update_each_iter(true, false);

// std::cout << " here 4" << std::endl;
    // Initialize V
    if (!Rf_isNull(noise_in["V"])) {
        V = Rcpp::as< VectorXd > (noise_in["V"]);
        prevV = V;
    } else {
        // V=h at init.
        // to-fix (for normal noise case)
        // sample_uncond_V(); // sample_uncond_V();
        for (int i=0; i < noise_type.size(); i++) {
            int n = V_size / n_noise;
            if (noise_type[i] == "normal") V.segment(i*n, n) = h.segment(i*n, n);
        }
    }


if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu(bool rao_blackwell) {
    // use conditional E(W|V,Y)
    VectorXd WW = (rao_blackwell) ? cond_W : W;
    VectorXd grad = VectorXd::Zero(n_theta_mu);
    if (n_theta_mu == 0 || fix_flag[latent_fix_theta_mu]) return grad;
    if ((n_noise == 1 && noise_type[0] == "normal") ||
        (n_noise == 2 && noise_type[0] == "normal" && noise_type[1] == "normal"))
        return grad;

    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

    // compute gradient g with V
    for (int l=0; l < n_theta_mu; l++) {
        grad(l) = (V-h).cwiseProduct(B_mu.col(l).cwiseQuotient(SV)).dot(getK()*WW - mu.cwiseProduct(V-h));
// std::cout << "grad mu  = " << grad(l) << std::endl;
// std::cout << "grad mu prior = " << PriorUtil::d_log_dens(prior_mu_type, prior_mu_param, theta_mu(l)) << std::endl;
        // add prior
        grad(l) += PriorUtil::d_log_dens(prior_mu_type, prior_mu_param, theta_mu(l));
    }

    return -grad;
}

// return the gradient wrt. theta, theta=log(sigma)
inline VectorXd Latent::grad_theta_sigma(bool rao_blackwell) {
    VectorXd WW = (rao_blackwell) ? cond_W : W;

    VectorXd grad = VectorXd::Zero(n_theta_sigma);
    if (fix_flag[latent_fix_theta_sigma]) return grad;

    // std::cout << "W = " << W << std::endl;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

    // (KW-mu(V-h)) / (sigma^2 V)
    VectorXd tmp = (getK()*WW - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseQuotient(SV);

    // grad = Bi(tmp * sigma ^ -2 - 1)
    grad = B_sigma.transpose() * (tmp - VectorXd::Ones(V_size));

    // add prior
    for (int l=0; l < n_theta_sigma; l++) {
        grad(l) += PriorUtil::d_log_dens(prior_sigma_type, prior_sigma_param, theta_sigma(l));
    }

    if (rao_blackwell) grad += rb_trace_sigma;
    return -grad;
}

inline VectorXd Latent::grad_theta_sigma_normal(bool rao_blackwell) {
    VectorXd WW = (rao_blackwell) ? cond_W : W;
    VectorXd V = VectorXd::Ones(V_size);
    VectorXd grad = VectorXd::Zero(n_theta_sigma_normal);

    // tmp = (KW - mu(V-h))^2 / V
    VectorXd tmp = (getK()*W).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
    // grad = Bi(tmp * sigma_normal ^ -2 - 1)
    VectorXd tmp1 = tmp.cwiseProduct(sigma_normal.array().pow(-2).matrix()) - VectorXd::Ones(V_size);
    grad += B_sigma_normal.transpose() * tmp1;

    // add prior
    for (int l=0; l < n_theta_sigma; l++) {
        grad(l) += PriorUtil::d_log_dens(prior_sigma_type, prior_sigma_param, theta_sigma(l));
    }
    // return - 1.0 / V_size * grad;
    return -grad;
}

VectorXd Latent::grad_theta_nu() {
// std::cout << "here nu" << std::endl;
    VectorXd grad = VectorXd::Zero(n_theta_nu);
    if (n_theta_nu == 0 || fix_flag[latent_fix_theta_nu]) return grad;

    if (noise_type[0] == "normal") return grad;
    if (n_noise == 1) {
        // single noise
        grad = NoiseUtil::grad_theta_nu(noise_type[0], B_nu, nu, V, prevV, h, single_V);
    } else { // same as single noise case
        // bv noise (does not allow for non-stationary nu)
        // 2 NIG case
        // n_theta_nu == 2 (1 for each field)
        grad = NoiseUtil::grad_theta_nu(noise_type[0], B_nu, nu, V, prevV, h, single_V);

        // int n = V_size / 2;
        // grad.head(1) = NoiseUtil::grad_theta_nu(noise_type[0], B_nu.block(0, 0, n, B_nu.cols()), theta_nu.head(1), V.segment(0, n), prevV.segment(0, n), h.segment(0, n), single_V);

        // if (share_V)
        //     grad(1) = grad(0);
        // else
        //     grad.tail(1) = NoiseUtil::grad_theta_nu(noise_type[1], B_nu.block(n, 0, n, B_nu.cols()), theta_nu.tail(1), V.segment(n, n), prevV.segment(n, n), h.segment(n, n), single_V);
    }

// std::cout << "here2" << std::endl;
    // add prior
//     for (int l=0; l < n_theta_nu; l++) {
// // std::cout << "grad  = " << grad(l) << std::endl;
// // std::cout << "grad of prior = " << PriorUtil::d_log_dens(prior_nu_type, prior_nu_param, nu(l)) << std::endl;
//         grad(l) -= PriorUtil::d_log_dens(prior_nu_type, prior_nu_param, nu(l));
//     }

// std::cout << "g_nu = " << grad << std::endl;
    return grad;
}

VectorXd Latent::grad_theta_K(bool rao_blackwell) {

    VectorXd grad = VectorXd::Zero(n_theta_K);
    if (fix_flag[latent_fix_theta_K]) return grad;

    VectorXd WW = (rao_blackwell) ? cond_W : W;
    if (numer_grad) {
        double val = logd_W_given_V(WW, getK(), mu, sigma, V);
        for (int i=0; i < n_theta_K; i++) {
            VectorXd tmp = theta_K;
            tmp(i) += eps;
            ope_addeps->update_K(tmp);
            double val_add_eps = logd_W_given_V(WW, ope_addeps->getK(), mu, sigma, V);
            grad(i) = (val_add_eps - val) / eps;
        }
    } else {
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
        VectorXd tmp = getK() * WW - mu.cwiseProduct(V-h);
        for (int j=0; j < n_theta_K; j++) {
            grad(j) = trace[j] - (ope->get_dK()[j] * WW).cwiseProduct(SV.cwiseInverse()).dot(tmp);
        }
    }
    // add trace term using RB
    if (rao_blackwell) grad -= rb_trace_K;

    // add prior
    for (int l=0; l < n_theta_K; l++) {
        grad(l) += PriorUtil::d_log_dens(prior_K_type, prior_K_param, theta_K(l));
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
        Rcpp::Named("theta_nu")     = theta_nu,
        Rcpp::Named("V")            = V,
        Rcpp::Named("W")            = W
    );
}

const VectorXd Latent::get_parameter() {
if (debug) std::cout << "Start latent get parameter"<< std::endl;
if (debug) std::cout << "fix_theta_sigma = " << 
    fix_flag[latent_fix_theta_sigma] << std::endl;

    VectorXd parameter = VectorXd::Zero(n_params);

    if (!fix_flag[latent_fix_theta_K])
        parameter.segment(0, n_theta_K) = theta_K;
    if (!fix_flag[latent_fix_theta_mu])
        parameter.segment(n_theta_K, n_theta_mu) = theta_mu;
    if (!fix_flag[latent_fix_theta_sigma])
        parameter.segment(n_theta_K+n_theta_mu, n_theta_sigma) = theta_sigma;
    if (!fix_flag[latent_fix_theta_nu])
        parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_theta_nu) = theta_nu;

    if (noise_type[0] == "normal_nig")
        parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_theta_nu, n_theta_sigma_normal) = theta_sigma_normal;

// if (debug) std::cout << "End latent get parameter"<< std::endl;
if (debug) {
    if (std::isnan(parameter(0))||std::isnan(-parameter(0)))
        throw std::runtime_error("isnan");
    std::cout << "parameter= " << parameter << std::endl;
}

    return parameter;
}

const VectorXd Latent::get_grad(bool rao_blackwell) {
if (debug) std::cout << "Start latent get grad"<< std::endl;
    VectorXd grad = VectorXd::Zero(n_params);

    // compute gradient of each parameter
    if (!fix_flag[latent_fix_theta_K])
        grad.segment(0, n_theta_K) = grad_theta_K(rao_blackwell);
    if (!fix_flag[latent_fix_theta_mu])
        grad.segment(n_theta_K, n_theta_mu) = grad_theta_mu(rao_blackwell);
    if (!fix_flag[latent_fix_theta_sigma])
        grad.segment(n_theta_K+n_theta_mu, n_theta_sigma) = grad_theta_sigma(rao_blackwell);
    if (!fix_flag[latent_fix_theta_nu])
        grad.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_theta_nu)  = grad_theta_nu();

    if (noise_type[0] == "normal_nig")
        grad.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_theta_nu, n_theta_sigma_normal) = grad_theta_sigma_normal(rao_blackwell);

if (debug) std::cout << "finish latent gradient"<< std::endl;
    return grad;
}

void Latent::set_parameter(const VectorXd& theta, bool update_dK) {
// if (debug) std::cout << "Start latent set parameter"<< std::endl;

    // nig, gal and normal+nig
    if (!fix_flag[latent_fix_theta_K])
        theta_K = theta.segment(0, n_theta_K);
    if (!fix_flag[latent_fix_theta_mu])
        theta_mu = theta.segment(n_theta_K, n_theta_mu);
    if (!fix_flag[latent_fix_theta_sigma])
        theta_sigma = theta.segment(n_theta_K+n_theta_mu, n_theta_sigma);
    if (!fix_flag[latent_fix_theta_nu])
        theta_nu = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma, n_theta_nu);
    
    if (noise_type[0] == "normal_nig") {
        theta_sigma_normal = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_theta_nu, n_theta_sigma_normal);
    }

    update_each_iter(false, update_dK);
if (debug) std::cout << "Finish latent set parameter"<< std::endl;
}

void Latent::sample_cond_V() {
    if (n_theta_nu == 0 || fix_flag[latent_fix_V]) return;
    prevV = V;

    // update b_inc (p,a_inc already built)
    b_inc = (getK() * W + mu.cwiseProduct(h)).cwiseQuotient(sigma).array().pow(2);

    int n = V_size / n_noise; // n is equal to h_size (of each mesh)
    // sample conditional V
    if (single_V) {
        if (share_V) {
            // type-G1 model (n_nu == 2, but keep same nu)
            double v1 = rGIG_cpp(-0.5-n, nu[0]+a_inc.dot(h),nu[0]+b_inc.dot(h.cwiseInverse()), latent_rng());
            V = v1 * h;
        } else {
            // type-G2 (n_nu==2) and also univariate single noise (n_nu==1)
            for (int i=0; i < n_noise; i++) {
                if (noise_type[i] == "normal") continue;

                double v1 = rGIG_cpp(
                    -0.5 - n*1.0/n_noise,
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
            for (int i=0; i < n_noise; i++) {
                if (noise_type[i] == "normal") continue;

                // sample conditional V
                V.segment(i*n, n) = rGIG_cpp((p_vec+p_inc).segment(i*n, n), (a_vec+a_inc).segment(i*n, n), (b_vec+b_inc).segment(i*n, n), latent_rng());
            }
        }

    }
}

void Latent::sample_uncond_V() {
    if (n_theta_nu == 0 || fix_flag[latent_fix_V]) return;
    prevV = V;
    int n = V_size / n_noise;

    // same logic as in simulation.R
    if (single_V) {
        for (int i=0; i < n_noise; i++) {
            double v;
            if (noise_type[i] == "nig" || noise_type[i] == "normal_nig")
                v = rGIG_cpp(-0.5, nu[i], nu[i], latent_rng());
            else if (noise_type[i] == "gal")
                v = rGIG_cpp(nu[i], 2*nu[i], 1e-14, latent_rng());
            V.segment(i*n, n) = v * h.segment(i*n, n);
        }
    } else {
        for (int i=0; i < n_noise; i++) {
            // sample unconditional V
            V.segment(i*n, n) = rGIG_cpp(p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), latent_rng());
        }
    }

    // keep both V same for bivariate noise with share_V case
    if (n_theta_nu == 2 && share_V) {
        V.segment(n, n) = V.segment(0, n);
    }
}

// at init, and after each set parameter
void Latent::update_each_iter(bool initialization, bool update_dK) {
// std::chrono::steady_clock::time_point startTime, endTime; startTime = std::chrono::steady_clock::now();
    if (!initialization) ope->update_K(theta_K);
// endTime = std::chrono::steady_clock::now(); std::cout << "K_update time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << std::endl;
// std::cout << " here 5" << std::endl;
// std::cout << "B_mu size = " << B_mu.rows() << " * " << B_mu.cols() << std::endl;
// std::cout << "t_mu size = " << theta_mu.size() << std::endl;
    mu = B_mu * theta_mu;
    sigma = (B_sigma * theta_sigma).array().exp();
    nu = (B_nu * theta_nu).array().exp();
    nu = nu.cwiseMax(nu_lower_bound);

    if (noise_type[0] == "normal_nig")
        sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

if (debug) std::cout << "update_each_iter" << std::endl;

    // update p,a,b, depend on nu, h
    if (n_theta_nu > 0) {
        for (int i=0; i < n_noise; i++) {
            if (noise_type[i] == "normal") continue;

            // update p_vec, a_vec, b_vec
            NoiseUtil::update_gig(noise_type[i], nu, p_vec, a_vec, b_vec, h, single_V);
            V = rGIG_cpp(p_vec, a_vec, b_vec, latent_rng());

            // update p_inc, a_inc
            p_inc = VectorXd::Constant(V_size, -0.5 * dim);
            a_inc = mu.cwiseQuotient(sigma).array().pow(2);
        }
    }

    // Update K and dK
    if (!fix_flag[latent_fix_theta_K]) {
        if (!numer_grad) {
            // Compute trace[i] = tr(K^-1 dK[i])
            ope->update_dK(theta_K);
            if (!zero_trace) {
                if (!symmetricK) {
                    if (W_size > 10) {
                        lu_solver_K.compute_LU(getK());
                        // lu_solver_K.compute_KTK(getK());
                        for (int i=0; i < n_theta_K; i++){
                            trace[i] = lu_solver_K.trace_num(ope->get_dK()[i]);
                        }
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
                    for (int i=0; i < n_theta_K; i++) {
                        trace[i] = chol_solver_K.trace_num(ope->get_dK()[i]);
                    }
                }
            }
        } else {
            // Update K numerical
            for (int i=0; i < n_theta_K; i++) {
                VectorXd tmp = theta_K;
                tmp(i) += eps;
                ope_addeps->update_K(tmp);
                if (update_dK) num_dK[i] = (ope_addeps->getK() - ope->getK()) / eps;
            }
        }
    }
if (debug) std::cout << "finish update_each_iter" << std::endl;
}

// for compute hessian
double Latent::log_density(const VectorXd& parameter, bool precond_K) {
    double logd_W = 0;

    VectorXd theta_K, theta_mu, theta_sigma, tmp_mu, tmp_sigma;
    if (precond_K) {
        // 1. pi(W|V)
        if (!fix_flag[latent_fix_theta_K])
            theta_K = parameter.head(n_theta_K);

        if (!fix_flag[latent_fix_theta_mu]) {
            theta_mu = parameter.segment(n_theta_K, n_theta_mu);
            tmp_mu = B_mu * theta_mu;
        } else {
            tmp_mu = mu;
        }

        if (!fix_flag[latent_fix_theta_sigma]) {
            theta_sigma = parameter.segment(n_theta_K + n_theta_mu, n_theta_sigma);
            tmp_sigma = (B_sigma * theta_sigma).array().exp();
        } else {
            tmp_sigma = sigma;
        }

        ope_precond->update_K(theta_K);
        SparseMatrix<double> K = ope_precond->getK();
        logd_W = logd_W_given_V(W, K, tmp_mu, tmp_sigma, prevV);
    } else {
        if (!fix_flag[latent_fix_theta_mu]) {
            theta_mu = parameter.head(n_theta_mu);
            tmp_mu = B_mu * theta_mu;
        } else {
            tmp_mu = mu;
        }

        if (!fix_flag[latent_fix_theta_sigma]) {
            theta_sigma = parameter.segment(n_theta_mu, n_theta_sigma);
            tmp_sigma = (B_sigma * theta_sigma).array().exp();
        } else {
            tmp_sigma = sigma;
        }
            
        logd_W = logd_KW_given_V(tmp_mu, tmp_sigma, prevV);
    }

    return -logd_W;
}

// log density of W|V
double Latent::logd_W_given_V(const VectorXd& W, const SparseMatrix<double>& K, const VectorXd& mu, const VectorXd& sigma, const VectorXd& V) {
    double l = 0;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

    VectorXd tmp = K * W - mu.cwiseProduct(V-h);
    if (K.rows() < 5) {
        MatrixXd Kd = K.toDense();
        l = log(Kd.diagonal().prod()) - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    } else {
        if (!symmetricK) {
            SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;

            // solver_Q.analyze(Q);
            solver_Q.compute(Q);
            l = 0.5 * solver_Q.logdet() - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        
            // SparseLU<SparseMatrix<double>> solver_Q_tmp;
            // solver_Q_tmp.analyzePattern(Q); 
            // solver_Q_tmp.factorize(Q);
            // l = 0.5 * solver_Q_tmp.logAbsDeterminant() - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        } else {
            chol_solver_K.compute(K);
// cholesky_solver solver_K_tmp;
// solver_K_tmp.analyze(K); solver_K_tmp.compute(K);
            l = chol_solver_K.logdet() - 0.5 * SV.array().log().sum()
                - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        }
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
MatrixXd Latent::precond(bool precond_K, double eps) {
    double numerical_eps = eps;

    VectorXd v;
    if (!precond_K) {
        v = VectorXd::Zero(n_theta_mu + n_theta_sigma);
        // v << theta_mu, theta_sigma;
        if (!fix_flag[latent_fix_theta_mu]) v.head(n_theta_mu) = theta_mu;
        if (!fix_flag[latent_fix_theta_sigma]) v.tail(n_theta_sigma) = theta_sigma;
    } else {
        v = VectorXd::Zero(n_params - n_theta_nu);
        // v << theta_K, theta_mu, theta_sigma;
        if (!fix_flag[latent_fix_theta_K]) v.head(n_theta_K) = theta_K;
        if (!fix_flag[latent_fix_theta_mu]) v.segment(n_theta_K, n_theta_mu) = theta_mu;
        if (!fix_flag[latent_fix_theta_sigma]) v.segment(n_theta_K+n_theta_mu, n_theta_sigma) = theta_sigma;
    }

    int n = v.size();
    MatrixXd num_hess_no_nu = VectorXd::Constant(n, 1.0).asDiagonal();

// which version of hessian to compute
if (!improve_hessian) {
// Forward difference (error is O(eps)
    // compute f_v = log_density(v + numerical_eps * e_i)
    double original_val = log_density(v, precond_K);
    VectorXd f_v (n);
    for (int i=0; i < n; i++) {
        VectorXd tmp_v = v; tmp_v(i) += numerical_eps;
        f_v(i) = log_density(tmp_v, precond_K);
    }
    
    // compute H_ij = d2 f / dxi dxj
    for (int i=0; i < n; i++) {
        for (int j=0; j <= i; j++) {
            VectorXd tmp_vij = v; tmp_vij(i) += numerical_eps; tmp_vij(j) += numerical_eps;
            double f_vij = log_density(tmp_vij, precond_K);
            num_hess_no_nu(i, j) = (f_vij - f_v(i) - f_v(j) + original_val) / (numerical_eps * numerical_eps);
        }
    }
} else {
// Central difference (error is O(eps^2)
    // compute fwd_fwd(i, j) = log_density(v + eps * e_i + eps * e_j)
	MatrixXd fwd_fwd (n, n);
	for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            VectorXd tmp = v; 
            tmp(i) += numerical_eps; tmp(j) += numerical_eps;
            fwd_fwd(i, j) = log_density(tmp, precond_K);
        }
	}

    MatrixXd fwd_bwd (n, n);
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            VectorXd tmp = v; 
            tmp(i) += numerical_eps; tmp(j) -= numerical_eps;
            fwd_bwd(i, j) = log_density(tmp, precond_K);
        }
    }       

    MatrixXd bwd_fwd (n, n);
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            VectorXd tmp = v; 
            tmp(i) -= numerical_eps; tmp(j) += numerical_eps;
            bwd_fwd(i, j) = log_density(tmp, precond_K);
        }
    }

    MatrixXd bwd_bwd (n, n);
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            VectorXd tmp = v; 
            tmp(i) -= numerical_eps; tmp(j) -= numerical_eps;
            bwd_bwd(i, j) = log_density(tmp, precond_K);
        }
    }

	// compute H_ij = d2 f / dxi dxj
	for (int i=0; i < n; i++) {
		for (int j=0; j <= i; j++) {
			num_hess_no_nu(i, j) = (fwd_fwd(i, j) - fwd_bwd(i, j) - bwd_fwd(i, j) + bwd_bwd(i, j)) / (4 * numerical_eps * numerical_eps);
		}
	}
}

	// fill in the lower triangular part
	for (int i=0; i < n; i++) {
		for (int j=0; j < i; j++) {
			num_hess_no_nu(j, i) = num_hess_no_nu(i, j);
		}
	}

    MatrixXd precond_full = MatrixXd::Zero(n_params, n_params);
    if (precond_K) {
        precond_full.topLeftCorner(n_params - n_theta_nu, n_params - n_theta_nu) = num_hess_no_nu;
    } else {
        // fill K part as diag(1)
        precond_full.topLeftCorner(n_theta_K, n_theta_K) =
            VectorXd::Constant(n_theta_K, V_size).asDiagonal();
        // fill the rest
        precond_full.block(n_theta_K, n_theta_K, n_theta_mu + n_theta_sigma, n_theta_mu + n_theta_sigma) = num_hess_no_nu;
    }

    // compute numerical hessian for theta_nu
    // watch-out! for some reason, hessian for theta_nu is getting too small
    if (fix_flag[latent_fix_theta_nu]) {
        precond_full.bottomRightCorner(n_theta_nu, n_theta_nu) = V_size * MatrixXd::Identity(n_theta_nu, n_theta_nu);
    } else {
        precond_full.bottomRightCorner(n_theta_nu, n_theta_nu) = V_size * MatrixXd::Identity(n_theta_nu, n_theta_nu);
        // precond_full.bottomRightCorner(n_theta_nu, n_theta_nu) = NoiseUtil::precond(numerical_eps, noise_type[0], prevV, h, B_nu, theta_nu, single_V);
    }

// std::cout << "print hessian: " << std::endl << precond_full << std::endl;
    return precond_full;
}

const SparseMatrix<double, 0, int>& Latent::get_dK(int i) {
    if (!numer_grad)
        return ope->get_dK()[i];
    else
        return num_dK[i];
}


// return log(pi(KW|V)) + log(pi(V))
double Latent::logd_W_V() {
    double logd_V = 0;
    if (n_noise == 1) {
        logd_V = NoiseUtil::log_density(noise_type[0], V, h, B_nu, theta_nu, single_V);
    } else if (n_noise == 2) {
        double logd_V1 = NoiseUtil::log_density(noise_type[0], V.head(V_size/2), h.head(V_size/2), B_nu.block(0, 0, V_size/2, B_nu.cols()), theta_nu, single_V);
        double logd_V2 = NoiseUtil::log_density(noise_type[1], V.tail(V_size/2), h.tail(V_size/2), B_nu.block(V_size/2, 0, V_size/2, B_nu.cols()), theta_nu, single_V);
        logd_V = logd_V1 + logd_V2;
    }

    return logd_W_given_V(W, getK(), mu, sigma, V) + logd_V;
}
