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

// handle specially, always return 0, but update Sigma
VectorXd Randeff::grad_theta_K() {
    VectorXd W = Ws[0];
    double V = vars[0].getV()(0);

    MatrixXd invSigma = Sigma.inverse();
    MatrixXd iSkroniS = kroneckerProduct(invSigma, invSigma);
    VectorXd UUt = vec(W * W.transpose());
    VectorXd dSigma_vech = 0.5 * Dd.transpose() * iSkroniS * (1/V * UUt - vec(Sigma));
    MatrixXd ddSigma = 0.5 * Dd.transpose() * iSkroniS * Dd;
    dSigma_vech = ddSigma.ldlt().solve(dSigma_vech);

    // test postive definite
    bool pos_def = false; double step = 1;
    while (!pos_def && step > 1e-3) {
        VectorXd vec_newSigma = Dd * (vech(Sigma) - dSigma_vech);
        MatrixXd newSigma = veci(vec_newSigma, W_size, W_size);
        SelfAdjointEigenSolver<MatrixXd> es(newSigma);
        if (es.eigenvalues().minCoeff() > 0) {
            pos_def = true;
            Sigma = newSigma;
        } else {
            step = step / 2;
    std::cout << "step = " << step << std::endl;
            dSigma_vech = step * dSigma_vech;
        }
    }
std::cout << "Sigma = " << Sigma << std::endl;

    return VectorXd::Zero(n_theta_K);
}

// called after set parameter
void Randeff::update_each_iter() {
    theta_K = vech(Sigma);
    K = getK(theta_K);
    for (int i=0; i < n_rep; i++)
        setSparseBlock(&K_rep, i*V_size, i*V_size, K);
}

void Randeff::sample_cond_V() {
    std::cout << "sample_cond_V in randeff" << std::endl;
    double a_inc_vec = mu.dot(Sigma.llt().solve(mu));

    for (int i=0; i < n_rep; i++) {
        VectorXd W = Ws[i];
        double b_inc_vec = (W+mu).dot(Sigma.llt().solve(W+mu));
        vars[i].sample_cond_V(a_inc_vec, b_inc_vec, W_size);
    }

    std::cout << "var.getV = " << vars[0].getV() << std::endl;
}