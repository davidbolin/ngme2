/*
    Matern2D model with stationary kappa:
        alpha is the smoothness parameter
        parameter_K(0) = kappa
        parameter_K(1) = kappa2
        K1 = kappa^2 * C + G
        K2 = kappa2^2 * C + G
        K = (D_2 *kron* I_n) * diag(c_1*K1,c_2*K2), where D_2 is dependence matrix
*/

/*
MATLAB version in multiK.m
    K0 = sparse(size(K2,1),size(K2,2));
    L = [c(1)*K1 K0; K0 c(2)*K2];
    K = kron(D,speye(size(K2,1)))*L;
    Ci = [Ci K0; K0 Ci];
*/
#include "../latent.h"
//model_list = latent_in
Matern2D::Matern2D(Rcpp::List& model_list, unsigned long seed)
: Latent(model_list, seed),
    G           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["G"])),
    C           (Rcpp::as< SparseMatrix<double,0,int> > (model_list["C"])),
    alpha       (Rcpp::as<int> (model_list["alpha"])),//assuming common alpha
    Cdiag       (C.diagonal())
    //TODO create new parameters in the latent model_list output
    theta       (Rcpp::as<double>(model_list["theta"])),
    rho         (Rcpp::as<double>(model_list["rho"])),
    c1          (Rcpp::as<double>(model_list["c1"])),
    c2          (Rcpp::as<double>(model_list["c2"])),
    typeG       (Rcpp::as<int>(model_list["typeG"]))//4 types for NIG SPDE model bivariate
{
std::cout << "begin Constructor of Matern " << std::endl;
    symmetricK = true;//TODO what does it control?
    parameter_K(0) = exp(parameter_K(0)); //exp(log(kappa1))
    parameter_K(1) = exp(parameter_K(1)); //exp(log(kappa2))
    theta_sigma(0) = exp(theta_sigma(0)); //exp(log(sigma1))
    theta_sigma(1) = exp(theta_sigma(1)); //exp(log(sigma2))

//setting rho and theta parameters accroding to the NIG model
    if(typeG == -1){
        rho = 0;
        theta = 0;
    }else if (typeG == 1)
    {
        rho = rho;
        theta = 0;
    }else if (typeG == 2)
    {
        rho = rho;
        theta = atan(rho);   
    }else { //typeG ==0, the most general case
        rho = rho;
        theta = theta;   
    }
    
    // Init K1 and K2
    K1 = getK(parameter_K(0));//= kappa1^2 * C + G
    K2 = getK(parameter_K(1));
    
    //initialize K to [G G;G G]
    SparseMatrix<double,0,int> K;
    //TODO is W_size = d??

    d = 2*G.rows();
    K.resize(d,d);
    setSparseBlock(&K,0,0, G);
    setSparseBlock(&K,d/2,d/2, G);
    setSparseBlock(&K,0,d/2, G);
    setSparseBlock(&K,d/2,0, G);

    // Init D dependence matrix 
    MatrixXd D(2,2);
    D(0,0) = cos(theta) + rho*sin(theta);
    D(0,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
    D(1,0) = sin(theta) - rho*cos(theta);
    D(1,1) = cos(theta)*pow(1+pow(rho,2),0.5);

    SparseMatrix<double,0,int> B;
    K1 = c1 * K1;//= c1 * (kappa1^2 * C + G)
    K2 = c2 * K2;
    //Set K from parameters by filling the four blocks:
    B = D(0,0)*K1;
    setSparseBlock_update(&K,0,0, B);
    B = D(1,1)*K2;
    setSparseBlock_update(&K,d/2,d/2, B);
    B = D(0,1)*K2;
    setSparseBlock_update(&K,0,d/2, B);
    B = D(1,0)*K1;
    setSparseBlock_update(&K,d/2,0, B);
    
    SparseMatrix<double> Q = K.transpose() * K;
    //TODO check if d= W_size?
    if (!use_iter_solver) {
        chol_solver_K.init(d, 0,0,0);
        chol_solver_K.analyze(K);
    } else {
        CG_solver_K.init(d,d,d, 0.5);//TODO why 0.5?
        CG_solver_K.analyze(K);
    }

    //Qsolver will only operate on blocks
    compute_trace();
    solver_Q.init(d, 0,0,0);
    solver_Q.analyze(Q); 

std::cout << "finish Constructor of Matern " << std::endl;
}

SparseMatrix<double> Matern2D::getK(VectorXd parameter_K) const {
    double kappa = parameter_K(0);

    int W_size = G.rows(); //TODO what is size of G, where does it come from?

    SparseMatrix<double> K_a (W_size, W_size);
        // VectorXd k2C = (kappa * kappa * Cdiag);
        // SparseMatrix<double> KCK = k2C.asDiagonal();
    SparseMatrix<double> KCK = kappa * kappa * C;

    // VectorXd kappas = VectorXd::Constant(W_size, parameter_K(0));
    // SparseMatrix<double> KCK (W_size, W_size);
    //     KCK = kappas.cwiseProduct(kappas).cwiseProduct(Cdiag).asDiagonal();

    if (alpha==2) {
        // K_a = T (G + KCK) C^(-1/2) -> Actually, K_a = C^{-1/2} (G+KCK), since Q = K^T K.
        K_a = (G + KCK);
    } else if (alpha==4) {
        // K_a = T (G + KCK) C^(-1) (G+KCK) C^(-1/2) -> Actually, K_a = C^{-1/2} (G + KCK) C^(-1) (G+KCK), since Q = K^T K.
        K_a = (G + KCK) *
            Cdiag.cwiseInverse().asDiagonal() * (G + KCK);
    } else {
        throw("alpha not equal to 2 or 4 is not implemented");
    }

    return K_a;
}

// stationary
SparseMatrix<double> Matern2D::get_dK(int index, VectorXd parameter_K) const {
//index, i, indicates the parameter_K[i]
//for dK wrt kappa1 index = 1; for kappa2 index =2;
    double kappa1 = parameter_K(0);
    double kappa2 = parameter_K(1);
    int W_size =  G.rows();
    SparseMatrix<double> dK1, dK2 (W_size, W_size);
    SparseMatrix<double> dK (2*W_size, 2*W_size);

    if (alpha==2)
        dK1 = 2*pow(kappa1,2)*C;
        dK2 = 2*pow(kappa2,2)*C;
    else if (alpha==4)
        dK1 = 4*pow(kappa1,2) * G + 4* pow(kappa1, 4) * C;
        dK2 = 4*pow(kappa2,2) * G + 4* pow(kappa2, 4) * C;
    else
        throw("alpha != 2 or 4");
    
    if (index == 1)
    {  
        double dc1 = -(alpha - 1)*c1;
        dK
    } else
    { //dL wrt kappa2
        double dc2 = -(alpha - 1)*c2;
        dK
    }
    //TODO can i initialize with zeros? without the matrix G
    B = D(0,0)*K1;
    setSparseBlock_update(&K,0,0, B);
    B = D(1,1)*K2;
    setSparseBlock_update(&K,d/2,d/2, B);
    B = D(0,1)*K2;
    setSparseBlock_update(&K,0,d/2, B);
    B = D(1,0)*K1;
    setSparseBlock_update(&K,d/2,0, B);
    return dK;
}

// compute numerical dK
void Matern::update_num_dK() {
    double kappa = parameter_K(0);
    double eps = 0.01;
    SparseMatrix<double> K_add_eps = pow(kappa + eps, 2) * C + G;
    dK = (K_add_eps - K) / eps;
}

VectorXd Matern::get_unbound_theta_K() const {
    assert (parameter_K.size() == 1);

    double th = k2th(parameter_K(0));
    return VectorXd::Constant(1, th);
}

// return length 1 vectorxd : grad_kappa * dkappa/dtheta
VectorXd Matern::grad_theta_K() {
    SparseMatrix<double> dK = get_dK_by_index(0);
    VectorXd V = getV();
    VectorXd SV = getSV();

    VectorXd kappa = parameter_K;
    double th = k2th(kappa(0));

    double da  = exp(th);
    double d2a = exp(th);

    double ret = 0;
    if (numer_grad) {
        // 1. numerical gradient
        ret = numerical_grad()(0);
    } else {
        // 2. analytical gradient and numerical hessian
        double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V).cwiseProduct(mu));
        double grad = trace - tmp;

        if (!use_precond) {
            ret = - grad * da / W_size;
        } else {
            VectorXd prevV = getPrevV();
            // compute numerical hessian
            SparseMatrix<double> K2 = getK_by_eps(0, eps);
            SparseMatrix<double> dK2 = get_dK_by_eps(0, 0, eps);

            // grad(x+eps) - grad(x) / eps
            VectorXd prevSV = getPrevSV();
            double grad2_eps = trace_eps - (dK2*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K2 * prevW +  (h - prevV).cwiseProduct(mu));
            double grad_eps  = trace - (dK*prevW).cwiseProduct(prevSV.cwiseInverse()).dot(K * prevW +  (h - prevV).cwiseProduct(mu));

            double hess = (grad2_eps - grad_eps) / eps;

            ret = (grad * da) / (hess * da * da + grad_eps * d2a);
        }
    }

    return VectorXd::Constant(1, ret);
}

void Matern::set_unbound_theta_K(VectorXd theta) {
    double kappa = th2k(theta(0));

    // update theta_K, K and dK
    parameter_K = VectorXd::Constant(1, kappa);
    K = getK(parameter_K);
    dK = get_dK(0, parameter_K);
    if (use_num_dK) {
        update_num_dK();
    }

    if (!numer_grad) compute_trace();
}

