// grad_theta_mu
// grad_theta_sigma
// function_K

#include "latent.h"
#include "operator.h"

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
        fix_flag[latent_fix_V]  = Rcpp::as<bool> (noise_in["fix_V"]);
        fix_flag[latent_fix_nu] = Rcpp::as<bool> (noise_in["fix_nu"]);
        if (!Rf_isNull(noise_in["V"])) {
            V = Rcpp::as< VectorXd > (noise_in["V"]);
            prevV = V;
        }

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

    update_each_iter(true);
if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
    VectorXd grad = VectorXd::Zero(n_theta_mu);
    if (fix_flag[latent_fix_theta_mu]) return grad;
    if (n_nu == 1 && noise_type[0] == "normal") return grad;

    double hess = 0;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

    hess += -(prevV-h).cwiseQuotient(prevSV).dot(prevV-h);
    for (int l=0; l < n_theta_mu; l++) {
        grad(l) += (V-h).cwiseProduct(B_mu.col(l).cwiseQuotient(SV)).dot(getK()*W - mu.cwiseProduct(V-h));
    }

    // not use hessian
    // if (V_size < 10)
    return - grad / sqrt(W_size);
    // else
        // return grad / (hess * n_rep);

// if (debug) std::cout << "KW" << K*W << std::endl;
// if (debug) std::cout << "V-h" << prevV-h << std::endl;
// if (debug) std::cout << "SV = " << getSV() << std::endl;

// if (debug) {
// std::cout << "grad of mu=" << grad <<std::endl;
// std::cout << "hess of mu=" << hess <<std::endl;
// }
    // return - grad / V_size;
}

// return the gradient wrt. theta, theta=log(sigma)
inline VectorXd Latent::grad_theta_sigma() {
    VectorXd grad = VectorXd::Zero(n_theta_sigma);
    if (fix_flag[latent_fix_theta_sigma]) return grad;

    // std::cout << "W = " << W << std::endl;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

    // tmp = (KW - mu(V-h))^2 / V
    VectorXd tmp = (getK()*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
    VectorXd tmp2 = tmp.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Ones(V_size);

    // grad = Bi(tmp * sigma ^ -2 - 1)
    grad = B_sigma.transpose() * tmp2;

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

const VectorXd Latent::get_parameter() const {
// if (debug) std::cout << "Start latent get parameter"<< std::endl;
    VectorXd parameter (n_params);

    parameter.segment(0, n_theta_K)                         = theta_K;
    parameter.segment(n_theta_K, n_theta_mu)                = theta_mu;
    parameter.segment(n_theta_K+n_theta_mu, n_theta_sigma)  = theta_sigma;
    parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_nu) = nu.array().log();
    if (noise_type[0] == "normal_nig")
        parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_nu, n_theta_sigma_normal) = theta_sigma_normal;
// if (debug) std::cout << "parameter= " << parameter << std::endl;
// if (debug) std::cout << "End latent get parameter"<< std::endl;
    return parameter;
}

const VectorXd Latent::get_grad() {
    VectorXd grad = VectorXd::Zero(n_params);

    grad.segment(0, n_theta_K) = grad_theta_K();
    grad.segment(n_theta_K, n_theta_mu) = grad_theta_mu();
    grad.segment(n_theta_K+n_theta_mu, n_theta_sigma) = grad_theta_sigma();
    grad.segment(n_theta_K+n_theta_mu+n_theta_sigma,n_nu)  = grad_theta_nu();
    if (noise_type[0] == "normal_nig")
        grad.segment(n_theta_K+n_theta_mu+n_theta_sigma+n_nu, n_theta_sigma_normal) = grad_theta_sigma_normal();

// std::cout << "g th k = " << grad_theta_K().transpose() << std::endl;
// std::cout << "g th mu = " << grad_theta_mu().transpose() << std::endl;
// std::cout << "g th sigma " << grad_theta_sigma().transpose() << std::endl;
// std::cout << "g th nu " << grad_theta_nu().transpose() << std::endl;

// DEBUG: checking grads
if (debug) {
    // std::cout << "gradient= " << grad << std::endl;
    // std::cout << "one latent gradient time " << since(grad1).count() << std::endl;
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

void Latent::update_each_iter(bool init) {
    // update mu and sigma
    mu = B_mu * theta_mu;
    sigma = (B_sigma * theta_sigma).array().exp();

    // update 1 by 1
    int n = V_size / n_nu;
    for (int i=0; i < n_nu; i++) {
        if (noise_type[i] == "gal")  {
            p_vec.segment(i*n, n) = h.segment(i*n, n) * nu(i);
            a_vec.segment(i*n, n) = VectorXd::Constant(n, nu(i) * 2);
            b_vec.segment(i*n, n) = VectorXd::Constant(n, 1e-14);
        } else { // nig or normal_nig
            p_vec.segment(i*n, n) = VectorXd::Constant(n, -0.5);
            a_vec.segment(i*n, n) = VectorXd::Constant(n, nu(i));
            b_vec.segment(i*n, n) = a_vec.segment(i*n, n).cwiseProduct(h.segment(i*n, n).cwiseProduct(h.segment(i*n, n)));
        }
        if (noise_type[i] == "normal") mu.segment(i*n, n).setZero();
    }

    if (noise_type[0] == "normal_nig")
        sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

    // update on K, dK, trace, bigK for sampling...
    if (!init) ope->update_K(theta_K);
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
            // update trace_eps if using hessian
            // if (use_precond) {
            //     SparseMatrix<double> dK_eps = get_dK_by_eps(i, 0, eps);
            //     if (!symmetricK) {
            //         lu_solver_K.computeKTK(K);
            //         trace[i] = lu_solver_K.trace(dK[i]);
            //     } else {
            //         chol_solver_K.compute(K);
            //         trace[i] = chol_solver_K.trace(dK[i]);
            //     }
            // }
        }
    }
}

double Latent::grad_theta_nu(
    const string& noise_type,
    double nu,
    const VectorXd& V,
    const VectorXd& prevV,
    const VectorXd& h
) const {
    if (fix_flag[latent_fix_nu] || noise_type == "normal") return 0;
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

void Latent::sample_V() {
    if (fix_flag[latent_fix_V]) return;
    prevV = V;

    if (n_nu == 1 && noise_type[0] != "normal")
        V = rGIG_cpp(p_vec, a_vec, b_vec, latent_rng());

    if (n_nu > 1) {
        int n = V_size / 2;
        if (noise_type[0] != "normal")
            V.segment(0, n) = rGIG_cpp(p_vec.segment(0, n), a_vec.segment(0, n), b_vec.segment(0, n), latent_rng());
        if (noise_type[1] != "normal")
            V.segment(n, n) = rGIG_cpp(p_vec.segment(n, n), a_vec.segment(n, n), b_vec.segment(n, n), latent_rng());
    }
}

void Latent::sample_cond_V() {
    if (fix_flag[latent_fix_V]) return;

    VectorXd a_inc_vec = mu.cwiseQuotient(sigma).array().pow(2);
    VectorXd b_inc_vec = (getK() * W + mu.cwiseProduct(h)).cwiseQuotient(sigma).array().pow(2);

    double dim = 1;
    VectorXd p_vec_new = p_vec - VectorXd::Constant(V_size, 0.5 * dim);
    VectorXd a_vec_new = a_vec + a_inc_vec;
    VectorXd b_vec_new = b_vec + b_inc_vec;

    prevV = V;
    // sample the part where is not normal
    if (n_nu == 1 && noise_type[0] != "normal")
        V = rGIG_cpp(p_vec_new, a_vec_new, b_vec_new, latent_rng());

    if (n_nu > 1) {
        int n = V_size / 2;
        if (noise_type[0] != "normal")
            V.segment(0, n) = rGIG_cpp(p_vec_new.segment(0, n), a_vec_new.segment(0, n), b_vec_new.segment(0, n), latent_rng());
        if (noise_type[1] != "normal")
            V.segment(n, n) = rGIG_cpp(p_vec_new.segment(n, n), a_vec_new.segment(n, n), b_vec_new.segment(n, n), latent_rng());
    }
}