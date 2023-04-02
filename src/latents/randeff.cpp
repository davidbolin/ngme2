// Notice here:
// for R interface, theta_K is directly the parameter_K of the Opeartor object
// for C interface, theta_K is f(parameter_K) such that it is unbounded
#include "../latent.h"

/*
    Random effect model:
        parameter_K = invSigma
*/

// W_size = V_size
// get_K_params, grad_K_params, set_K_params, output
Randeff::Randeff(const Rcpp::List& model_list, unsigned long seed)
    : Latent(model_list, seed),
    Sigma      (Rcpp::as<MatrixXd> (model_list["Sigma"])),
    Dd         (duplicatematrix(W_size))
{
if (debug) std::cout << "Begin Constructor of Randeff1" << std::endl;
    // theta_K = vech(Sigma)

    // Init K and Q
    K = getK(theta_K);
if (debug) std::cout << "Begin Constructor of Randeff2" << std::endl;
    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);
if (debug) std::cout << "K = " << K << std::endl;

    // Init Q
    SparseMatrix<double> Q = K.transpose() * K;
    solver_Q.init(W_size, 0,0,0);
    solver_Q.analyze(Q);
if (debug) std::cout << "End Constructor of Randeff1" << std::endl;
}

//     VectorXd getM() const {
//         double V = var.getV()(0);
// V = 1;
//         return V*(V-1) * Sigma * mu;
//     }

// K = sqrt(invSigma)
SparseMatrix<double> Randeff::getK(const VectorXd& Sigma_vech) const {
    // std::cout << "Sigma_vech = " << Sigma_vech << std::endl;
    VectorXd Sigma_vec = Dd * Sigma_vech;
    // std::cout << "Sigma_vec = " << Sigma_vec << std::endl;
    MatrixXd Sigma = veci(Sigma_vec, W_size, W_size);
    SparseMatrix<double> K = Sigma.inverse().llt().matrixL().toDenseMatrix().sparseView();
    return K;
}

SparseMatrix<double> Randeff::get_dK(int index, const VectorXd& alpha) const {
    throw std::runtime_error("Not implemented yet");
}

VectorXd Randeff::grad_theta_K() {
    VectorXd W = Ws[0];

    MatrixXd invSigma = Sigma.inverse();
    MatrixXd iSkroniS = kroneckerProduct(invSigma, invSigma);
    VectorXd UUt = vec(W * W.transpose());
    VectorXd dSigma_vech = 0.5 * Dd.transpose() * iSkroniS  * (UUt - vec(Sigma));
std::cout << "dSigma_vech = " << dSigma_vech << std::endl;
    return dSigma_vech;
}

// called after set parameter
void Randeff::update_each_iter() {
    K = getK(theta_K);

    // update Sigma
    VectorXd Sigma_vec = Dd * theta_K;
    Sigma = veci(Sigma_vec, W_size, W_size);

    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);
}

void Randeff::sample_cond_V() {
    VectorXd a_inc_vec = mu * Sigma.llt().solve(mu);

    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        VectorXd b_inc_vec = (W+mu) * Sigma.llt().solve(W+mu);
        vars[i].sample_cond_V(a_inc_vec, b_inc_vec, W_size);
    }
}