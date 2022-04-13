#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;

class AR : public Latent {
    
public:
    AR(Rcpp::List ar1_in) 
    : Latent(ar1_in)
    {
        
        // read init kappa
        Rcpp::List init_value = Rcpp::as<Rcpp::List> (ar1_in["init_value"]);
        double kappa = Rcpp::as<double>  (init_value["kappa"]);

        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        ope = new GC(ope_in, kappa);
        
        // Init K
        solver_K.init(n_reg, 0,0,0);
        solver_K.analyze(getK());
        compute_trace();
    }

    void sample_cond_V() {
        var->sample_cond_V(getK(), W, h, mu, sigma);
    };
    

    double th2a(double th) const {
        return (-1 + 2*exp(th) / (1+exp(th)));
    }
    double a2th(double a) const {
        return (log((-1-a)/(-1+a)));
    }

    double get_theta_kappa() const {
        double a = getKappa();
        return a2th(a);
    }

    void set_theta_kappa(double v) {
        double a = th2a(v);
        setKappa(a);
    }

    // grad_kappa * dkappa/dtheta
    double grad_theta_kappa() {
        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
        VectorXd V = getV();
        VectorXd prevV = getPrevV();
        
        VectorXd SV = pow(sigma, 2) * V;
        VectorXd prevSV = pow(sigma, 2) * prevV;
        
        double a = getKappa();
        double th = a2th(a);

        double da  = 2 * (exp(th) / pow(1+exp(th), 2));
        double d2a = -2 * (exp(th) * (-1+exp(th)) / pow(1+exp(th), 3));

        double tmp = (dK*W).cwiseProduct(SV.cwiseInverse()).dot(K * W + (h - V) * mu);
        double grad = trace - tmp;

        tmp = (dK*W).cwiseProduct(prevSV.cwiseInverse()).dot(dK * W);
        double hess = -trace2 - tmp;

std::cout << "******* grad of kappa: " << grad << std::endl;   
std::cout << "******* hess of kappa: " << hess << std::endl;   

        return (grad * da) / (hess * da * da + grad * d2a) * n_reg;
    }

};
