// grad_theta_mu
// grad_theta_sigma
// function_K

#include "latent.h"
#include "operator.h"

// K is V_size * W_size matrix
Latent::Latent(const Rcpp::List& model_list, unsigned long seed) :
    latent_rng    (seed),
    model_type    (Rcpp::as<string>     (model_list["model"])),
    noise_type    (Rcpp::as<string>     (model_list["noise_type"])),
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
    nu            (1)
    // var           (Var(Rcpp::as<Rcpp::List> (model_list["noise"]), latent_rng()))
{
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
        if (noise_type=="normal_nig") B_sigma_normal  = Rcpp::as< MatrixXd >    (noise_in["B_sigma_normal"]);

        n_theta_mu    =   (B_mu.cols());
        n_theta_sigma =   (B_sigma.cols());
        if (noise_type=="normal_nig") n_theta_sigma_normal =   (B_sigma_normal.cols());

        theta_mu = Rcpp::as< VectorXd >    (noise_in["theta_mu"]);
        theta_sigma = Rcpp::as< VectorXd > (noise_in["theta_sigma"]);
        if (noise_type=="normal_nig") theta_sigma_normal = Rcpp::as< VectorXd > (noise_in["theta_sigma_normal"]);

        nu = Rcpp::as<double> (noise_in["nu"]);
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

    // Init noise
    if (noise_type == "normal") {
        fix_flag[latent_fix_theta_mu] = 1; // no mu need
        fix_flag[latent_fix_nu] = 1;
    }

    // init Q
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
    update_each_iter();
if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
    VectorXd grad = VectorXd::Zero(n_theta_mu);
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
    if (numer_grad) {
        double val = function_K(getK());
        for (int i=0; i < n_theta_K; i++) {
            VectorXd tmp = theta_K;
            tmp(i) += eps;
            ope_add_eps->update_K(tmp);
            SparseMatrix<double> K_add_eps = ope_add_eps->getK();
            double val_add_eps = function_K(K_add_eps);
            grad(i) = - (val_add_eps - val) / eps;
// std::cout << "num_g = " << grad(i)  << std::endl;
        }
    } else {
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

        VectorXd tmp = getK() * W - mu.cwiseProduct(V-h);
        for (int j=0; j < n_theta_K; j++) {
            grad(j) = trace[j] - (ope->get_dK()[j] * W).cwiseProduct(SV.cwiseInverse()).dot(tmp);
            grad(j) = - grad(j) / W_size;
        }
    }

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

    if (noise_type == "normal") {
        parameter.segment(0, n_theta_K)              = theta_K;
        parameter.segment(n_theta_K, n_theta_sigma)  = theta_sigma;
    } else {
    // nig, gal, and nig+normal
        parameter.segment(0, n_theta_K)                         = theta_K;
        parameter.segment(n_theta_K, n_theta_mu)                = theta_mu;
        parameter.segment(n_theta_K+n_theta_mu, n_theta_sigma)  = theta_sigma;
        parameter(n_theta_K+n_theta_mu+n_theta_sigma)           = log(nu);
    if (noise_type == "normal_nig")
        parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma+1, n_theta_sigma_normal) = theta_sigma_normal;
    }

// if (debug) std::cout << "parameter= " << parameter << std::endl;
// if (debug) std::cout << "End latent get parameter"<< std::endl;
    return parameter;
}

const VectorXd Latent::get_grad() {
    VectorXd grad = VectorXd::Zero(n_params);

    if (noise_type == "normal") {
        if (!fix_flag[latent_fix_theta_K])
            grad.segment(0, n_theta_K) = grad_theta_K();
        if (!fix_flag[latent_fix_theta_sigma])
            grad.segment(n_theta_K, n_theta_sigma) = grad_theta_sigma();
    } else {
        if (!fix_flag[latent_fix_theta_K])
            grad.segment(0, n_theta_K) = grad_theta_K();
        if (!fix_flag[latent_fix_theta_mu])
            grad.segment(n_theta_K, n_theta_mu) = grad_theta_mu();
        if (!fix_flag[latent_fix_theta_sigma])
            grad.segment(n_theta_K+n_theta_mu, n_theta_sigma) = grad_theta_sigma();
        grad(n_theta_K+n_theta_mu+n_theta_sigma)  = grad_theta_nu();

        if (noise_type == "normal_nig")
            grad.segment(n_theta_K+n_theta_mu+n_theta_sigma+1, n_theta_sigma_normal) = grad_theta_sigma_normal();
    }

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
    if (noise_type == "normal") {
        theta_K = theta.segment(0, n_theta_K);
        theta_sigma = theta.segment(n_theta_K, n_theta_sigma);
        sigma = (B_sigma * theta_sigma).array().exp();
    } else {
        // nig, gal and normal+nig
        theta_K = theta.segment(0, n_theta_K);
        theta_mu = theta.segment(n_theta_K, n_theta_mu);
        theta_sigma = theta.segment(n_theta_K+n_theta_mu, n_theta_sigma);
        double log_nu = (theta(n_theta_K+n_theta_mu+n_theta_sigma));
        nu = exp(log_nu);

        if (noise_type == "normal_nig") {
            theta_sigma_normal = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma+1, n_theta_sigma_normal);
        }
    }
    update_each_iter();
}

void Latent::update_each_iter() {
    int n = V_size;
    if (noise_type == "gal")  {
        p_vec = h * nu;
        a_vec = VectorXd::Constant(n, nu * 2);
        b_vec = VectorXd::Constant(n, 1e-14);
    } else { // nig or normal_nig
        p_vec = VectorXd::Constant(n, -0.5);
        a_vec = VectorXd::Constant(n, nu);
        b_vec = a_vec.cwiseProduct(h.cwiseProduct(h));
    }

    mu = B_mu * theta_mu;
    sigma = (B_sigma * theta_sigma).array().exp();
    if (noise_type=="normal_nig")
        sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

    // update on K, dK, trace, bigK for sampling...
    ope->update_K(theta_K);
    if (!numer_grad && W_size == V_size) {
        ope->update_dK(theta_K);
        if (!zero_trace) {
            for (int i=0; i < n_theta_K; i++) {
                if (!symmetricK) {
                    lu_solver_K.computeKTK(getK());
                    trace[i] = lu_solver_K.trace(ope->get_dK()[i]);
                } else {
                    chol_solver_K.compute(getK());
                    trace[i] = chol_solver_K.trace(ope->get_dK()[i]);
                }
    // std::cout << "trace = " << trace[i] << std::endl;
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
}

double Latent::grad_theta_nu() {
    if (fix_flag[latent_fix_nu] || noise_type == "normal") return 0;
    int n = V_size;
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

void Latent::sample_V() {
    if (fix_flag[latent_fix_V] || noise_type == "normal") return;

    prevV = V;
    V = rGIG_cpp(p_vec, a_vec, b_vec, latent_rng());
}

void Latent::sample_cond_V() {
    if (fix_flag[latent_fix_V] || noise_type == "normal") return;

    VectorXd a_inc_vec = mu.cwiseQuotient(sigma).array().pow(2);
    VectorXd b_inc_vec = (getK() * W + mu.cwiseProduct(h)).cwiseQuotient(sigma).array().pow(2);

    double dim = 1;
    VectorXd p_vec_new = p_vec - VectorXd::Constant(V_size, 0.5 * dim);
    VectorXd a_vec_new = a_vec + a_inc_vec;
    VectorXd b_vec_new = b_vec + b_inc_vec;

    prevV = V;
    V = rGIG_cpp(p_vec_new, a_vec_new, b_vec_new, latent_rng());
}