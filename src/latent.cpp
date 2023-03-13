// grad_theta_mu
// grad_theta_sigma
// function_K
// numerical_grad

#include "latent.h"

// K is V_size * W_size matrix
Latent::Latent(const Rcpp::List& model_list, unsigned long seed) :
    latent_rng    (seed),
    model_type    (Rcpp::as<string>     (model_list["model"])),
    noise_type    (Rcpp::as<string>     (model_list["noise_type"])),
    debug         (Rcpp::as<bool>       (model_list["debug"])),

    n_rep         (Rcpp::as<int>        (model_list["n_rep"])),
    W_size        (Rcpp::as<int>        (model_list["W_size"])),
    V_size        (Rcpp::as<int>        (model_list["V_size"])),
    n_params      (Rcpp::as<int>        (model_list["n_params"])),
    theta_K       (Rcpp::as<VectorXd>   (model_list["theta_K"])),
    n_theta_K     (Rcpp::as<int>        (model_list["n_theta_K"])),

    K            (V_size, W_size),
    dK           (V_size, W_size),
    K_rep        (n_rep * V_size, n_rep * W_size),

    trace         (0),
    trace_eps     (0),
    eps           (0.01),

    // W             (W_size),
    // prevW         (W_size),
    h             (Rcpp::as< VectorXd >                     (model_list["h"])), //same length as V_size
    A             (Rcpp::as< SparseMatrix<double,0,int> >   (model_list["A"])),

    Ws            (n_rep),
    prevWs        (n_rep),
    vars          (n_rep),
    // var           (Var(Rcpp::as<Rcpp::List> (model_list["noise"]), latent_rng())),

    theta_K_traj  (theta_K.size()),

    left (nullptr),
    right (nullptr)
{
if (debug) std::cout << "Begin constructor of latent" << std::endl;

    // init Vs : create V for each replicate
    for (int i=0; i < n_rep; i++) {
        vars[i] = Var(Rcpp::as<Rcpp::List> (model_list["noise"]), latent_rng());
    }

    // setting the seed
    // latent_rng.seed(seed);

    // read from ngme.model
    fix_flag[latent_fix_theta_K] = Rcpp::as<bool>    (model_list["fix_theta_K"]);

    // read the control variable
    Rcpp::List control_f = Rcpp::as<Rcpp::List> (model_list["control"]);
        use_precond     = Rcpp::as<bool>        (control_f["use_precond"] );
        numer_grad      = Rcpp::as<bool>        (control_f["numer_grad"]) ;
        eps             = Rcpp::as<double>      (control_f["eps"]) ;
        use_iter_solver = Rcpp::as<bool>        (control_f["use_iter_solver"]);

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

        mu = (B_mu * theta_mu);
        sigma = (B_sigma * theta_sigma).array().exp();
        if (noise_type=="normal_nig") sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();

    const int n_nu = 1;

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

    theta_mu_traj.resize(n_theta_mu);
    theta_sigma_traj.resize(n_theta_sigma);
    if (noise_type=="normal_nig") theta_sigma_normal_traj.resize(n_theta_sigma_normal);
    record_traj();

if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
// if (debug) std::cout << "Start mu gradient"<< std::endl;
    VectorXd grad = VectorXd::Zero(n_theta_mu);
    double hess = 0;

    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
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
    return grad / (hess * n_rep);

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
        VectorXd V = vars[i].getV();
        VectorXd prevV = vars[i].getPrevV();
        VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
        VectorXd prevSV = sigma.array().pow(2).matrix().cwiseProduct(prevV);

        // tmp = (KW - mu(V-h))^2 / V
        VectorXd tmp = (K*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
        VectorXd tmp1 = tmp.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Ones(V_size);

        tmp2 += tmp1;
// std::cout << "W size = " << W.size() <<std::endl;
// std::cout << "V size = " << V.size() <<std::endl;
// std::cout << "mu size = " << mu.size() <<std::endl;
// std::cout << "sigma size = " << sigma.size() <<std::endl;
// std::cout << "h size = " << h.size() <<std::endl;
// std::cout << "V_size = " << V_size <<std::endl;
// std::cout << "W_size = " << W_size <<std::endl;
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
    return l / n_rep;
}

// function_K(params += ( 0,0,eps,0,0) )
double Latent::function_K(VectorXd& theta_K) {
    SparseMatrix<double> K = getK(theta_K);

    return function_K(K);
}

// numerical gradient for K parameters
VectorXd Latent::numerical_grad() {
    SparseMatrix<double> K = getK(theta_K);
    VectorXd grad = VectorXd::Zero(n_theta_K);

    double val = function_K(K);
    // iterate every parameter
    for (int i=0; i < n_theta_K; i++) {
        SparseMatrix<double> K_add_eps = getK_by_eps(i, eps);
        double val_add_eps = function_K(K_add_eps);
        double num_g = (val_add_eps - val) / eps;

        if (!use_precond) {
            grad(i) = - num_g / W_size;
        } else {
            SparseMatrix<double> K_minus_eps = getK_by_eps(i, -eps);
            double val_minus_eps = function_K(K_minus_eps);
            double num_hess = (val_minus_eps + val_add_eps - 2*val) / pow(eps, 2);
            grad(i) = num_g / num_hess;
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

    // Rcpp::List trajecotry = Rcpp::List::create(
    //     Rcpp::Named("theta_K")            = theta_K_traj,
    //     Rcpp::Named("theta_mu")           = theta_mu_traj,
    //     Rcpp::Named("theta_sigma")        = theta_sigma_traj,
    //     Rcpp::Named("theta_sigma_normal")        = theta_sigma_normal_traj,
    //     Rcpp::Named("nu")            = nu_traj
    // );
    // out.attr("trajectory") = trajecotry;
    return out;
}


// // numerical gradient for K parameters
// VectorXd Latent::numerical_grad() {
// std::cout << "start numerical gradient" <<std::endl;
//     int n_theta_K =   get_n_params();
//     VectorXd params = get_parameter();
//     double val = function_K(params);

//     VectorXd grad (n_theta_K);
//     for (int i=0; i < n_theta_K; i++) {
//         VectorXd params_add_eps = params;
//             params_add_eps(i) += eps;
//         double val_add_eps = function_K(params_add_eps);
//         double num_g = (val_add_eps - val) / eps;

//         if (!use_precond) {
//             grad(i) = - num_g / W_size;
//             // grad(i) = - num_g;
//         } else {
//             VectorXd params_minus_eps = params;
//                 params_minus_eps(i) -= eps;
//             double val_minus_eps = function_K(params_minus_eps);
//             double num_hess = (val_minus_eps + val_add_eps - 2*val) / pow(eps, 2);
//             grad(i) = num_g / num_hess;
//         }
//     }
//     // return grad;
//     return grad;
// }