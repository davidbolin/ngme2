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
    ope           (OperatorFactory::create(Rcpp::as<Rcpp::List> (model_list["operator"]))),
    theta_K       (Rcpp::as<VectorXd>  (model_list["theta_K"])),
    n_theta_K     (theta_K.size()),
    symmetricK    (ope->is_symmetric()),

    n_rep         (1),
    W_size        (Rcpp::as<int>        (model_list["W_size"])),
    V_size        (Rcpp::as<int>        (model_list["V_size"])),
    n_params      (Rcpp::as<int>        (model_list["n_params"])),

    K            (V_size, W_size),
    dK           (n_theta_K),

    trace         (n_theta_K, 0.0),
    eps           (0.001),

    // W             (W_size),
    // prevW         (W_size),
    h             (ope->get_h()),
    A             (Rcpp::as< SparseMatrix<double,0,int> >   (model_list["A"])),

    Ws            (n_rep),
    prevWs        (n_rep),
    vars          (n_rep)
    // var           (Var(Rcpp::as<Rcpp::List> (model_list["noise"]), latent_rng())),
{
    debug = true;
std::cout << "Begin constructor of latent" << std::endl;

    // init Vs : create V for each replicate
    for (int i=0; i < n_rep; i++) {
        vars[i] = Var(Rcpp::as<Rcpp::List> (model_list["noise"]), latent_rng());
    }

    // read the control variable
    Rcpp::List control_f = Rcpp::as<Rcpp::List> (model_list["control"]);
        use_precond     = Rcpp::as<bool>        (control_f["use_precond"] );
        numer_grad      = Rcpp::as<bool>        (control_f["numer_grad"]) ;
        eps             = Rcpp::as<double>      (control_f["eps"]) ;

    // construct from ngme.noise
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

    // init W
    for (int i=0; i < n_rep; i++) {
        if (model_list["W"] != R_NilValue) {
            Ws[i] = Rcpp::as< VectorXd > (model_list["W"]);
        } else {
            Ws[i] = VectorXd::Zero(W_size);
        }
        prevWs[i] = Ws[i];
    }
    fix_flag[latent_fix_W] = Rcpp::as<bool> (model_list["fix_W"]); // fixW

    // About noise
    // if (fix_flag[latent_fix_V]) var.fixV();
    if (vars[0].get_noise_type() == "normal") {
        fix_flag[latent_fix_theta_mu] = 1; // no mu need
    }

    // init K and Q
    K = ope->getK(theta_K);

    if (!symmetricK) {
        lu_solver_K.init(W_size, 0,0,0);
        lu_solver_K.analyze(K);
    } else {
        chol_solver_K.init(W_size,0,0,0);
        chol_solver_K.analyze(K);
    }

    SparseMatrix<double> Q = K.transpose() * K;
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);

    update_each_iter();
if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
    VectorXd grad = VectorXd::Zero(n_theta_mu);
    double hess = 0;

    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
// if (debug) std::cout << "Kw = "<< K * W << std::endl;
        VectorXd V = vars[i].getV();
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
        VectorXd prevV = vars[i].getPrevV();
        VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

        hess += -(prevV-h).cwiseQuotient(prevSV).dot(prevV-h);
        for (int l=0; l < n_theta_mu; l++) {
            grad(l) += (V-h).cwiseProduct(B_mu.col(l).cwiseQuotient(SV)).dot(K*W - mu.cwiseProduct(V-h));
        }
    }
    hess /= n_rep;

    // not use hessian
    // if (V_size < 10)
        return - grad / (sqrt(W_size) * n_rep);
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
    VectorXd tmp2 = VectorXd::Zero(V_size);

    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
    // std::cout << "W = " << W << std::endl;
        VectorXd V = vars[i].getV();
        VectorXd prevV = vars[i].getPrevV();
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
        VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

        // tmp = (KW - mu(V-h))^2 / V
        VectorXd tmp = (K*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
        VectorXd tmp1 = tmp.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Ones(V_size);

        tmp2 += tmp1;
    }

    // grad = Bi(tmp * sigma ^ -2 - 1)
    grad = B_sigma.transpose() * (tmp2 / n_rep);

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
    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        VectorXd tmp = (K*W).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
        // grad = Bi(tmp * sigma_normal ^ -2 - 1)
        VectorXd tmp1 = tmp.cwiseProduct(sigma_normal.array().pow(-2).matrix()) - VectorXd::Ones(V_size);
        grad += B_sigma_normal.transpose() * tmp1;
    }
    return - 1.0 / V_size * grad;
}

double Latent::function_K(SparseMatrix<double>& K) {
    double l = 0;
    for (int i = 0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        VectorXd V = vars[i].getV();
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

        VectorXd tmp = K * W - mu.cwiseProduct(V-h);
        if (K.rows() < 5) {
            MatrixXd Kd = K.toDense();
            l += log(Kd.diagonal().prod()) - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
        } else {
            if (!symmetricK) {
                SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
                solver_Q.compute(Q);
                l += 0.5 * solver_Q.logdet()
                    - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
            } else {
                chol_solver_K.compute(K);
                l += chol_solver_K.logdet()
                    - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
            }
        }
    }
    return l / n_rep;
}

VectorXd Latent::grad_theta_K() {
    VectorXd grad = VectorXd::Zero(n_theta_K);
    if (numer_grad) {
        double val = function_K(K);
        for (int i=0; i < n_theta_K; i++) {
            VectorXd tmp = theta_K;
            tmp(i) += eps;
            SparseMatrix<double> K_add_eps = ope->getK(tmp);
            double val_add_eps = function_K(K_add_eps);
// std::cout << " theta_ K  =" << theta_K << std::endl;
// std::cout << " K = " << K << std::endl;
// std::cout << "Kadd eps = " << K_add_eps << std::endl;
            double num_g = (val_add_eps - val) / eps;
            grad(i) = - num_g / W_size;
        }
    } else {
        for (int i=0; i < n_rep; i++) {
            VectorXd W = Ws[i];
// std::cout << " W = " << W << std::endl;
            VectorXd V = vars[i].getV();
            VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);

            VectorXd tmp = K * W - mu.cwiseProduct(V-h);
            for (int j=0; j < n_theta_K; j++) {
                grad(j) = trace[j] - (dK[j] * W).cwiseProduct(SV.cwiseInverse()).dot(tmp);
                grad(j) = - grad(j) / W_size;
            }
        }
    }

    return grad;
}

Rcpp::List Latent::output() const {
    // compute mean of V and W
    VectorXd W = getW(); VectorXd V = getV();
    VectorXd meanV = VectorXd::Zero(V_size);
    VectorXd meanW = VectorXd::Zero(W_size);
    for (int i=0; i < n_rep; i++) {
        meanV += V.segment(i*V_size, V_size);
        meanW += W.segment(i*W_size, W_size);
    }
    meanV /= n_rep; meanW /= n_rep;

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("model")        = model_type,
        Rcpp::Named("noise_type")   = noise_type,
        Rcpp::Named("theta_K")      = theta_K, // same parameterization as input
        Rcpp::Named("theta_mu")     = theta_mu,
        Rcpp::Named("theta_sigma")  = theta_sigma,
        Rcpp::Named("theta_sigma_normal")  = theta_sigma_normal,
        Rcpp::Named("nu")           = vars[0].get_nu(),  // gives eta > 0, not log(eta)
        Rcpp::Named("V")            = meanV,
        Rcpp::Named("W")            = meanW
    );

    return out;
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
        parameter(n_theta_K+n_theta_mu+n_theta_sigma)           = vars[0].get_log_nu();
    if (noise_type == "normal_nig")
        parameter.segment(n_theta_K+n_theta_mu+n_theta_sigma+1, n_theta_sigma_normal) = theta_sigma_normal;
    }

// if (debug) std::cout << "parameter= " << parameter << std::endl;
// if (debug) std::cout << "End latent get parameter"<< std::endl;
    return parameter;
}

const VectorXd Latent::get_grad() {
// if (debug) std::cout << "Start latent gradient"<< std::endl;
// auto grad1 = std::chrono::steady_clock::now();
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
        for (int i=0; i < n_rep; i++) vars[i].set_log_nu(log_nu); // for each replicate

        if (noise_type == "normal_nig") {
            theta_sigma_normal = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma+1, n_theta_sigma_normal);
        }

        update_each_iter();
    }
}

void Latent::update_each_iter() {
        mu = (B_mu * theta_mu);
        sigma = (B_sigma * theta_sigma).array().exp();
        if (noise_type=="normal_nig") sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

        // update on K, dK, trace, bigK for sampling...
        K = ope->getK(theta_K);
        if (!numer_grad && W_size == V_size) {
            for (int i=0; i < n_theta_K; i++) {
                dK[i] = ope->get_dK(i, theta_K);
            }
            if (!zero_trace) {
                for (int i=0; i < n_theta_K; i++) {
                    if (!symmetricK) {
                        lu_solver_K.computeKTK(K);
                        trace[i] = lu_solver_K.trace(dK[i]);
                    } else {
                        chol_solver_K.compute(K);
                        trace[i] = chol_solver_K.trace(dK[i]);
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