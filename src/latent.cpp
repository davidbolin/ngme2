// grad_theta_mu
// grad_theta_sigma
// function_K
// function_kappa
// numerical_grad

#include "latent.h"

// K is V_size * W_size matrix
Latent::Latent(Rcpp::List& model_list, unsigned long seed) :
    latent_rng    (seed),
    model_type    (Rcpp::as<string>     (model_list["model"])),
    noise_type    (Rcpp::as<string>     (model_list["noise_type"])),
    debug         (Rcpp::as<bool>       (model_list["debug"])),
    W_size        (Rcpp::as<int>        (model_list["W_size"])),
    V_size        (Rcpp::as<int>        (model_list["V_size"])),
    parameter_K   (Rcpp::as<VectorXd>   (model_list["theta_K"])),
    n_theta_K     (Rcpp::as<int>        (model_list["n_theta_K"])),

    trace         (0),
    trace_eps     (0),
    eps           (0.001),

    W             (W_size),
    prevW         (W_size),
    h             (Rcpp::as< VectorXd >                     (model_list["h"])), //same length as V_size
    A             (Rcpp::as< SparseMatrix<double,0,int> >   (model_list["A"])),

    var           (Var(Rcpp::as<Rcpp::List> (model_list["noise"]), latent_rng())),

    theta_K_traj  (parameter_K.size())
{
if (debug) std::cout << "Begin constructor of latent" << std::endl;

    // setting the seed
    // latent_rng.seed(seed);

    // read from ngme.model
    fix_flag[latent_fix_theta_K] = Rcpp::as<bool>    (model_list["fix_theta_K"]);
    fix_flag[latent_fix_W]       = Rcpp::as<bool>    (model_list["fix_W"]);

    // read the control variable
    Rcpp::List control_f = Rcpp::as<Rcpp::List> (model_list["control"]);
        use_precond     = Rcpp::as<bool>        (control_f["use_precond"] );
        numer_grad      = Rcpp::as<bool>        (control_f["numer_grad"]) ;
        eps             = Rcpp::as<double>      (control_f["eps"]) ;
        use_iter_solver = Rcpp::as<bool>        (control_f["use_iter_solver"]);

    string noise_type = Rcpp::as<string>     (model_list["noise_type"]);

    // construct from ngme.noise
    Rcpp::List noise_in = Rcpp::as<Rcpp::List> (model_list["noise"]);
        fix_flag[latent_fix_theta_mu]     = Rcpp::as<bool>  (noise_in["fix_theta_mu"]);
        fix_flag[latent_fix_theta_sigma]  = Rcpp::as<bool>  (noise_in["fix_theta_sigma"]);
        fix_flag[latent_fix_theta_V]      = Rcpp::as<bool>  (noise_in["fix_theta_V"]);
        fix_flag[latent_fix_V]            = Rcpp::as<bool>  (noise_in["fix_V"]);

        B_mu     = Rcpp::as< MatrixXd >    (noise_in["B_mu"]);
        B_sigma  = Rcpp::as< MatrixXd >    (noise_in["B_sigma"]);
        n_theta_mu    =   (B_mu.cols());
        n_theta_sigma =   (B_sigma.cols());

        theta_mu_traj.resize(n_theta_mu);
        theta_sigma_traj.resize(n_theta_sigma);

        theta_mu = Rcpp::as< VectorXd >    (noise_in["theta_mu"]);
        set_theta_mu(theta_mu);
        theta_sigma = Rcpp::as< VectorXd > (noise_in["theta_sigma"]);
        set_theta_sigma(theta_sigma);

    const int n_theta_V = 1;
    n_params = n_theta_K + n_theta_mu + n_theta_sigma + n_theta_V;

    if (model_list["W"] != R_NilValue) {
        W = Rcpp::as< VectorXd >    (model_list["W"]);
        prevW = W;
    }

    // About noise
    if (fix_flag[latent_fix_V]) var.fixV();
    if (var.get_noise_type() == "normal") {
        fix_flag[latent_fix_theta_mu] = 1; // no mu need
    }

if (debug) std::cout << "End constructor of latent" << std::endl;
}

VectorXd Latent::grad_theta_mu() {
// if (debug) std::cout << "Start mu gradient"<< std::endl;
    VectorXd result(n_theta_mu);

    // VectorXd inv_V = V.cwiseInverse();
    // VectorXd prevV = getPrevV();
    // VectorXd prev_inv_V = prevV.cwiseInverse();
    VectorXd V = getV();

    VectorXd grad (n_theta_mu);
    for (int l=0; l < n_theta_mu; l++) {
        grad(l) = (V-h).cwiseProduct( B_mu.col(l).cwiseQuotient(getSV()) ).dot(K*W - mu.cwiseProduct(V-h));
    }

    result = - 1.0 / V_size * grad;

// if (debug) {
// std::cout << "grad of mu=" << grad <<std::endl;
// std::cout << "hess of mu=" << hess <<std::endl;
// }
    return result;
}

// return the gradient wrt. theta, theta=log(sigma)
inline VectorXd Latent::grad_theta_sigma() {
    VectorXd V = getV();
    VectorXd prevV = getPrevV();

    VectorXd result(n_theta_sigma);
    // double msq = (K*W - mu.cwiseProduct(V-h)).cwiseProduct(V.cwiseInverse()).dot(K*W - mu(0)*(V-h));
    // VectorXd vsq = (K*W - mu.cwiseProduct(V-h)).array().pow(2);
    VectorXd vsq = (K*W - mu.cwiseProduct(V-h)).array().pow(2).matrix().cwiseProduct(V.cwiseInverse());
    VectorXd grad (n_theta_sigma);
    // for (int l=0; l < n_theta_sigma; l++) {
    //     VectorXd tmp1 = vsq.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Constant(V_size, 1);
    //     VectorXd tmp2 = B_sigma.col(l).cwiseProduct(tmp1);
    //     grad(l) = tmp2.sum();
    // }

    // vector manner
    VectorXd tmp1 = vsq.cwiseProduct(sigma.array().pow(-2).matrix()) - VectorXd::Constant(V_size, 1);
    grad = B_sigma.transpose() * tmp1;

    VectorXd prev_vsq = (K*prevW - mu.cwiseProduct(prevV-h)).array().pow(2).matrix().cwiseProduct(prevV.cwiseInverse());
    MatrixXd hess (n_theta_sigma, n_theta_sigma);
    VectorXd tmp3 = -2*prev_vsq.cwiseProduct(sigma.array().pow(-2).matrix());

    hess = B_sigma.transpose() * tmp3.asDiagonal() * B_sigma;

    result = - 1.0 / V_size * grad;

    // result = hess.llt().solve(grad);

    return result;
}

double Latent::function_K(SparseMatrix<double> K) {
    VectorXd V = getV();
    VectorXd SV = getSV();

    SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
    VectorXd tmp = K * W - mu.cwiseProduct(V-h);

    double l;
    if (!symmetricK) {
        solver_Q.compute(Q);
        l = 0.5 * solver_Q.logdet()
               - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    } else {
        chol_solver_K.compute(K);
        l = chol_solver_K.logdet()
            - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    }
    return l;
}

// function_K(params += ( 0,0,eps,0,0) )
double Latent::function_K(VectorXd parameter_K) {
    SparseMatrix<double> K = getK(parameter_K);

    VectorXd V = getV();
    VectorXd SV = getSV();

    SparseMatrix<double> Q = K.transpose() * SV.cwiseInverse().asDiagonal() * K;
    VectorXd tmp = K * W - mu.cwiseProduct(V-h);

    double l;
    if (!symmetricK) {
        solver_Q.compute(Q);
        l = 0.5 * solver_Q.logdet()
               - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    } else {
        chol_solver_K.compute(K);
        l = chol_solver_K.logdet()
            - 0.5 * tmp.cwiseProduct(SV.cwiseInverse()).dot(tmp);
    }
    return l;
}

// numerical gradient for K parameters
VectorXd Latent::numerical_grad() {
std::cout << "start numerical gradient" <<std::endl;
    double val = function_K(K);

    VectorXd grad (n_theta_K);
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
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("model")        = model_type,
        Rcpp::Named("noise_type")   = noise_type,
        Rcpp::Named("theta_K")      = parameter_K, // same parameterization as input
        Rcpp::Named("theta_mu")     = theta_mu,
        Rcpp::Named("theta_sigma")  = theta_sigma,
        Rcpp::Named("theta_V")      = var.get_theta_V(),  // gives eta > 0, not log(eta)
        Rcpp::Named("V")            = getV(),
        Rcpp::Named("W")            = W
    );
    Rcpp::List trajecotry = Rcpp::List::create(
        Rcpp::Named("theta_K")            = theta_K_traj,
        Rcpp::Named("theta_mu")           = theta_mu_traj,
        Rcpp::Named("theta_sigma")        = theta_sigma_traj,
        Rcpp::Named("theta_V")            = theta_V_traj
    );
    out.attr("trajectory") = trajecotry;
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