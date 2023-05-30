#include "latent.h"

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

        nu = Rcpp::as<VectorXd> (noise_in["nu"]);
        n_nu = nu.size();

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
        sample_uncond_V(); sample_uncond_V();
    }

if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
    VectorXd grad = VectorXd::Zero(n_theta_mu);
    if (fix_flag[latent_fix_theta_mu]) return grad;
    if (n_nu == 1 && noise_type[0] == "normal") return grad;
    if (n_nu == 2 && noise_type[0] == "normal" && noise_type[1] == "normal") return grad;

    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

    for (int l=0; l < n_theta_mu; l++) {
        grad(l) += (V-h).cwiseProduct(B_mu.col(l).cwiseQuotient(SV)).dot(getK()*W - mu.cwiseProduct(V-h));
    }
    double hess = -(prevV-h).cwiseQuotient(prevSV).dot(prevV-h);

    if (V_size < 10)
        return - grad / sqrt(W_size);
    else
        return grad / hess;
        // return - grad / sqrt(V_size);
        // return - 1.0 / V_size * grad;
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

    return - 1.0 / V_size * grad;

    // for (int l=0; l < n_theta_sigma; l++) {
    //     VectorXd tmp1 = tmp.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Ones(V_size);
    //     VectorXd tmp2 = B_sigma.col(l).cwiseProduct(tmp1);
    //     grad(l) = tmp2.sum();
    // }

// compute hessian
    // VectorXd prevV = getPrevV();
    // VectorXd prev_tmp = (K*prevW - mu.cwiseProduct(prevV-h)).array().pow(2).matrix().cwiseProduct(prevV.cwiseInverse());
    // MatrixXd hess (n_theta_sigma, n_theta_sigma);
    // VectorXd tmp3 = -2*prev_tmp.cwiseProduct(sigma.array().pow(-2).matrix());
    // hess = B_sigma.transpose() * tmp3.asDiagonal() * B_sigma;

    // return hess.llt().solve(grad);
}

inline VectorXd Latent::grad_theta_sigma_normal() {
    VectorXd V = VectorXd::Ones(V_size);
    VectorXd grad = VectorXd::Zero(n_theta_sigma_normal);

    // tmp = (KW - mu(V-h))^2 / V
    VectorXd tmp = (getK()*W).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
    // grad = Bi(tmp * sigma_normal ^ -2 - 1)
    VectorXd tmp1 = tmp.cwiseProduct(sigma_normal.array().pow(-2).matrix()) - VectorXd::Ones(V_size);
    grad += B_sigma_normal.transpose() * tmp1;

    return - 1.0 / V_size * grad;
}

// pi(W|V)
double Latent::function_K(const SparseMatrix<double>& K) {
    double l = 0;
    VectorXd SV = getSV();

    VectorXd tmp = K * W - mu.cwiseProduct(V-h);
    if (K.rows() < 5) {
        MatrixXd Kd = K.toDense();
        l = log(Kd.diagonal().prod()) - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    } else {
        if (!symmetricK) {
            SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
            solver_Q.compute(Q);
            l = 0.5 * solver_Q.logdet() - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        } else {
            chol_solver_K.compute(K);
            l = chol_solver_K.logdet() - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        }
    }
    // normalize
// std::cout << " l = " << l / W_size << std::endl;
    return l / W_size;
}

VectorXd Latent::grad_theta_K() {
// std::cout << "K = " << K << std::endl;
    VectorXd grad = VectorXd::Zero(n_theta_K);
    if (fix_flag[latent_fix_theta_K]) return grad;

    if (numer_grad) {
        double val = function_K(getK());
        for (int i=0; i < n_theta_K; i++) {
            VectorXd tmp = theta_K;
            tmp(i) += eps;
            ope_add_eps->update_K(tmp);
            SparseMatrix<double> K_add_eps = ope_add_eps->getK();
            double val_add_eps = function_K(K_add_eps);
            grad(i) = - (val_add_eps - val) / eps;
        }
    } else {
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

        VectorXd tmp = getK() * W - mu.cwiseProduct(V-h);
        for (int j=0; j < n_theta_K; j++) {
            grad(j) = trace[j] - (ope->get_dK()[j] * W).cwiseProduct(SV.cwiseInverse()).dot(tmp);
            grad(j) = - grad(j) / W_size;
        }
    }
// std::cout << "grad_K = " << grad.transpose() << std::endl;
    return grad;
}

Rcpp::List Latent::output() const {
    return  Rcpp::List::create(
        Rcpp::Named("model")        = model_type,
        Rcpp::Named("noise_type")   = noise_type,
        Rcpp::Named("theta_K")      = theta_K, // same parameterization as input
        Rcpp::Named("theta_mu")     = theta_mu,
        Rcpp::Named("theta_sigma")  = theta_sigma,
        Rcpp::Named("theta_sigma_normal")  = theta_sigma_normal,
        Rcpp::Named("nu")           = nu,  // gives eta > 0, not log(eta)
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
    VectorXd grad = VectorXd::Zero(n_params);

    // compute gradient of each parameter
    grad.segment(0, n_theta_K) = grad_theta_K();
    grad.segment(n_theta_K, n_theta_mu) = grad_theta_mu();
    grad.segment(n_theta_K+n_theta_mu, n_theta_sigma) = grad_theta_sigma();
    grad.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_nu)  = grad_theta_nu();
    if (noise_type[0] == "normal_nig")
        grad.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_nu, n_theta_sigma_normal) = grad_theta_sigma_normal();
std::cout << "gradgg = " << grad << std::endl;
// DEBUG: checking grads
if (debug) {
    if (abs(grad.segment(0, n_theta_K).mean()) > 1)
    std::cout << "g_K is large, g_K = " << grad.segment(0, n_theta_K).mean() << std::endl;

    if (abs(grad.segment(n_theta_K, n_theta_mu).mean()) > 1)
    std::cout << "g_mu is large, g_mu = " << grad.segment(n_theta_K, n_theta_mu).mean() << std::endl;

    if (abs(grad.segment(n_theta_K+n_theta_mu, n_theta_sigma).mean()) > 1)
    std::cout << "g_sigma is large, g_sigma = " << grad.segment(n_theta_K+n_theta_mu, n_theta_sigma).mean() << std::endl;

    if (abs(grad.segment(n_theta_K+n_theta_mu+n_theta_sigma, n_nu).mean()) > 1)
    std::cout << "g_nu is large, g_nu = " << grad.segment(n_theta_K+n_theta_mu+n_theta_sigma, n_nu).mean() << std::endl;
}
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
    if (fix_flag[latent_fix_V]) return;
    prevV = V;

    // update b_inc (p,a_inc already built)
    b_inc = (getK() * W + mu.cwiseProduct(h)).cwiseQuotient(sigma).array().pow(2);

    int n = V_size / n_nu; // n is equal to h_size (of each mesh)
    // sample conditional V
    if (single_V) {
        if (share_V) {
            // type-G1 model (n_nu == 2, but keep same nu)
            double v = rGIG_cpp(-0.5-n, nu[0]+a_inc.dot(h),nu[0]+b_inc.dot(h.cwiseInverse()), latent_rng());
            V = v * VectorXd::Ones(V_size);
        } else {
            // type-G2 (n_nu==2) and also univariate single noise (n_nu==1)
            for (int i=0; i < n_nu; i++) {
                double v = rGIG_cpp(
                    -0.5 - n*1.0/n_nu,
                    nu[i] + a_inc.segment(i*n, n).dot(h.segment(i*n, n)),
                    nu[i] + b_inc.segment(i*n, n).dot(h.segment(i*n, n).cwiseInverse()),
                    latent_rng()
                );
                V.segment(i*n, n) = v * VectorXd::Ones(n);
            }
        }
    } else {
        if (share_V) {
            // type-G3 (n_nu == 2, but keep same nu)
            V = rGIG_cpp(VectorXd::Constant(V_size, -1.5), a_vec+a_inc, b_vec+b_inc, latent_rng());
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
    if (fix_flag[latent_fix_V]) return;
    prevV = V;

    int n = V_size / n_nu;
    if (single_V) {
        for (int i=0; i < n_nu; i++) {
            double v;
            if (noise_type[i] == "nig" || noise_type[i] == "normal_nig")
                v = rGIG_cpp(-0.5, nu[i], nu[i], latent_rng());
            else if (noise_type[i] == "gal")
                v = rGIG_cpp(-0.5, nu[i], nu[i], latent_rng());
            V.segment(i*n, n) = v * VectorXd::Ones(n);
        }
        if (n_nu == 2 && share_V) {
            V.segment(n, n) = V.segment(0, n);
        }
    } else {
        for (int i=0; i < n_nu; i++) {
            // sample unconditional V
            V.segment(i*n, n) = rGIG_cpp(p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), latent_rng());
        }
        if (n_nu == 2 && share_V) {
            V.segment(n, n) = V.segment(0, n);
        }
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
    int n = V_size / n_nu;
    for (int i=0; i < n_nu; i++) {
        if (noise_type[i] == "normal") continue;
        // update p_vec, a_vec, b_vec
        NoiseUtil::update_gig(noise_type[i], nu(i),
        p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), h.segment(i*n, n));
        V.segment(i*n, n) = rGIG_cpp(p_vec.segment(i*n, n), a_vec.segment(i*n, n), b_vec.segment(i*n, n), latent_rng());

        // update p_inc, a_inc
        p_inc = VectorXd::Constant(V_size, -0.5 * dim);
        a_inc = mu.cwiseQuotient(sigma).array().pow(2);
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