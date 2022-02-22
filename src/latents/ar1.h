#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"

class AR : public Latent {
    lu_sparse_solver solver_K;
    
public:
    AR(Rcpp::List ar1_in) 
    : Latent(ar1_in)
    {
        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        ope = new GC(ope_in);

        solver_K.init(n_reg, 0,0,0);
        solver_K.analyze(getK());
    }

    void sample_cond_V() {
        var->sample_cond_V(getK(), W, h, mu, sigma);
    };

    double _grad_kappa() {
        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
        VectorXd V = getV();
        solver_K.compute(K);

        double lhs = solver_K.trace(dK); // tr(dK * K^-1)
// std::cout << "factorized lhs=" << lhs << std::endl; 
        // 2. Compute the rest
        double rhs = W.transpose() * dK.transpose() * 
                    (VectorXd::Constant(n_reg, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * mu);
// std::cout << "factorized rhs=" << rhs << std::endl; 
        return (1.0/n_reg) * (rhs - lhs);
    }

};
