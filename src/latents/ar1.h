#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"

class AR : public Latent {
    SparseLU<SparseMatrix<double> > solver;
    // lu_sparse_solver K_solver;

    // unsigned n_obs;
    // VectorXd theta_K, theta_m, theta_V, Mu, W, grad;

    // SparseMatrix<double,0,int> A;
    // Operator K;
    // Var var;
public:
    AR(){}
    AR(Rcpp::List ar1_in) 
    {
        A                 = Rcpp::as< SparseMatrix<double,0,int> > (ar1_in["A"]);
        n_obs             = Rcpp::as< unsigned > (ar1_in["n"]);
        double theta      = Rcpp::as< double > (ar1_in["theta_init"]);
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        Rcpp::List var_in = Rcpp::as<Rcpp::List> (ar1_in["var_in"]);
        
        // Init
        Theta.resize(1); Grad.resize(1); Mu.resize(n_obs); W.resize(n_obs);
        Mu = VectorXd::Constant(n_obs, 0);

        // W comes from block level
        W = VectorXd::Constant(n_obs, 1);

        Theta << theta;

        K = GC(ope_in);
        init_var(var_in);

        // K_solver.init(n_obs, 1, 1, 0);
        solver.analyzePattern(getK());
        // compute_grad();
    }

    void compute_grad() {
        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
        SparseMatrix<double> d2K = get_d2K();
        VectorXd Mu = getMu();
        VectorXd W = getW();
        VectorXd V = getV();
        VectorXd h = VectorXd::Constant(V.size(), 1);

        // K_solver.analyze(K);
        // K_solver.compute(K);
        solver.factorize(K);
        
        // grad = K_solver.trace(dK) - W.transpose() * dK.transpose() * (VectorXd::Constant(V.size(), 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * Mu);
        
        // 1. To compute tr(dK * K^-1)
        SparseMatrix<double> I (n_obs, n_obs); 
        I.setIdentity();
        double g = (dK * solver.solve(I)).eval().diagonal().sum();
        
        // 2. Compute the rest // Debug
        double rhs = W.transpose() * dK.transpose() * (VectorXd::Constant(n_obs, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * Mu);
        Grad << g - rhs;


std::cout << "Mu=" << Mu;  std::cout << "W=" << W;  std::cout << "V=" << V;
std::cout << "K=" << K; std::cout << "dK=" << dK;  // std::cout << "W=" << W;  std::cout << "V=" << V;

// std::cout << g;
    }

};