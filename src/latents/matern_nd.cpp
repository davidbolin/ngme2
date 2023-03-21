#include "../latent.h"

class Matern_nd : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    VectorXd Cdiag;
public:
    Matern_nd(const Rcpp::List& model_list, unsigned long seed);
    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
    VectorXd grad_theta_K();
    void update_each_iter();
    void update_num_dK();
    double th2k(double th) const {return exp(th);}
    double k2th(double k) const {return log(k);}

    VectorXd unbound_to_bound_K(const VectorXd& theta_K) const {
        VectorXd kappa (1);
        kappa(0) = th2k(theta_K(0));
        return kappa;
    }
    VectorXd bound_to_unbound_K(const VectorXd& kappa) const {
        VectorXd theta_K (1);
        theta_K(0) = k2th(kappa(0));
        return theta_K;
    }

    VectorXd grad_theta_K(
        SparseMatrix<double>& K,
        SparseMatrix<double>& dK,
        vector<VectorXd>& Ws,
        vector<VectorXd>& prevWs,
        vector<Var>& vars,
        const VectorXd& mu,
        const VectorXd& sigma,
        const VectorXd& h,
        double trace,
        int W_size
    );
};

Matern_nd::Matern_nd(const Rcpp::List& model_list, unsigned long seed)
  : Latent(model_list, seed) {
  alpha = model_list["alpha"];
}

SparseMatrix<double> Matern_nd::getK(const VectorXd& alpha) const {
  // build Q
  SparseMatrix<double> Q;
  // build Dl
  SparseMatrix<double> Dl;

}
