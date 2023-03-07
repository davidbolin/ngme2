/*
    latent class:
    basically, providing get, grad, set (theta_K, ...)
        1. theta_K is unbounded
        2. theta_mu
        3. theta_sigma
        4. nu
*/

#ifndef NGME_LATANT_H
#define NGME_LATANT_H

// #include<Eigen/IterativeLinearSolvers>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <random>
#include <cmath>

#include "include/timer.h"
#include "include/solver.h"
#include "var.h"

using std::exp;
using std::log;
using std::pow;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

enum Latent_fix_flag {
    latent_fix_theta_K,
    latent_fix_W,
    latent_fix_theta_mu,
    latent_fix_theta_sigma
};
const int LATENT_FIX_FLAG_SIZE = 4;

class Latent {
protected:
    std::mt19937 latent_rng;
    string model_type, noise_type;
    bool debug;
    int n_rep, W_size, V_size, n_params, n_var {1}; // n_params=n_theta_K + n_theta_mu + n_theta_sigma + n_var

    // operator K related
    VectorXd theta_K;
    int n_theta_K;

    bool use_num_dK {false};
    SparseMatrix<double, 0, int> K, dK, d2K; // size = V_size * W_size
    SparseMatrix<double, 0, int> K_rep;  // n_rep of size of K (for sampling W)

    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0};

    bool use_precond {false}, numer_grad {false};
    bool symmetricK {false};

    // mu and sigma, and sigma_normal (special case when using nig_normal case)
    MatrixXd B_mu, B_sigma, B_sigma_normal;
    VectorXd theta_mu, theta_sigma, theta_sigma_normal;

    // mu = Bmu * theta_mu
    // sigma = exp(Bsigma * theta_sigma)
    VectorXd mu, sigma, sigma_normal;
    int n_theta_mu, n_theta_sigma, n_theta_sigma_normal;

    // eps for numerical gradient.
    double trace, trace_eps, eps;

    // VectorXd W, prevW, h;
    VectorXd h;
    SparseMatrix<double,0,int> A;

    vector<VectorXd> Ws, prevWs;
    vector<Var> vars;

    // Var var;

    // solver
    cholesky_solver chol_solver_K;
    lu_sparse_solver lu_solver_K;
    bool use_iter_solver {false};
    // iterative_solver CG_solver_K;

    cholesky_solver solver_Q; // Q = KT diag(1/SV) K

    // record trajectory
    vector<vector<double>> theta_K_traj;
    vector<vector<double>> theta_mu_traj;
    vector<vector<double>> theta_sigma_traj;
    vector<vector<double>> theta_sigma_normal_traj;
    vector<double>   nu_traj;
public:
    Latent(const Rcpp::List&, unsigned long seed);
    virtual ~Latent() {}

    // change of variable
    virtual VectorXd bound_to_unbound_K(const VectorXd&) const {return theta_K;}
    virtual VectorXd unbound_to_bound_K(const VectorXd&) const {return theta_K;}
    virtual void update_each_iter()=0;

    /*  1 Model itself   */
    int get_W_size() const                  {return n_rep * W_size; }
    int get_V_size() const                  {return n_rep * V_size; }
    int get_n_params() const                {return n_params; }
    SparseMatrix<double, 0, int>& getA()    {return A; }

    // concat W of different replicate
    const VectorXd getW() const {
        VectorXd W(n_rep * W_size);
        for (int i=0; i < n_rep; i++) {
            W.segment(i*W_size, W_size) = Ws[i];
        }
        return W;
    }

    const VectorXd getPrevW()  const {
        VectorXd W(n_rep * W_size);
        for (int i=0; i < n_rep; i++) {
            W.segment(i*W_size, W_size) = prevWs[i];
        }
        return W;
    }

    // newW is of size n_rep * W_size
    void setW(const VectorXd& newW) {
        if (!fix_flag[latent_fix_W]) {
            // prevW = W; W = newW;
            for (int i=0; i < n_rep; i++) {
                prevWs[i] = Ws[i];
                Ws[i] = newW.segment(i*W_size, W_size);
            }
        }
    }

    void setPrevW(const VectorXd& W) {
        for (int i=0; i < n_rep; i++) {
            prevWs[i] = W.segment(i*W_size, W_size);
        }
    }

    // return of dim V_size
    // const VectorXd getVMean() const {
    //     VectorXd meanV(V_size);
    //     for (int i=0; i < n_rep; i++) {
    //         meanV += vars[i].getV();
    //     }
    //     return meanV / n_rep;
    // }

    // used in block model
    VectorXd getMean() const {
        VectorXd Mean (n_rep * V_size);
        for (int i=0; i < n_rep; i++) {
            VectorXd V = vars[i].getV();
            Mean.segment(i*V_size, V_size) = mu.cwiseProduct(V - h);
        }
        return Mean;
    }

    /*  2 Variance component   */
    const VectorXd getV() const {
        VectorXd V(V_size * n_rep);
        for (int i=0; i < n_rep; i++) {
            V.segment(i*V_size, V_size) = vars[i].getV();
        }
        return V;
    }
    const VectorXd getPrevV() const {
        VectorXd V(V_size * n_rep);
        for (int i=0; i < n_rep; i++) {
            V.segment(i*V_size, V_size) = vars[i].getPrevV();
        }
        return V;
    }
    const VectorXd getSV() const {
        VectorXd SV(V_size * n_rep);
        for (int i=0; i < n_rep; i++) {
            SV.segment(i*V_size, V_size) = sigma.array().pow(2).matrix().cwiseProduct(vars[i].getV());
        }
        return SV;
    }
    void setPrevV(const VectorXd& V) {
        for (int i=0; i < n_rep; i++) {
            vars[i].setPrevV(V.segment(i*V_size, V_size));
        }
    }

    // VectorXd getSV() const {
    //     VectorXd V = getV();
    //     return (sigma.array().pow(2).matrix().cwiseProduct(V));
    // }
    // VectorXd getPrevSV() const {
    //     VectorXd prevV = getPrevV();
    //     return (sigma.array().pow(2).matrix().cwiseProduct(prevV));
    // }

    void sample_V() {
        for (int i=0; i < n_rep; i++) {
            vars[i].sample_V();
        }
    }

    void sample_cond_V() {
        VectorXd a_inc_vec = mu.cwiseQuotient(sigma).array().pow(2);
        for (int i=0; i < n_rep; i++) {
            VectorXd W = Ws[i];
            VectorXd tmp = (K * W + mu.cwiseProduct(h));
            VectorXd b_inc_vec = tmp.cwiseQuotient(sigma).array().pow(2);
            vars[i].sample_cond_V(a_inc_vec, b_inc_vec);
        }
    }

    /*  3 Operator component   */
    SparseMatrix<double, 0, int>& getK()  {
        return K_rep;
    }

    // get K/dK using different parameter
    virtual SparseMatrix<double, 0, int> getK(const VectorXd&) const=0;
    virtual SparseMatrix<double, 0, int> get_dK(int, const VectorXd&) const=0;

    // variants of getK and get_dK
    // get_dK wrt. paramter_K[i]
    SparseMatrix<double, 0, int> get_dK_by_index(int index) const {
        return get_dK(index, theta_K);
    }

    // param(pos) += eps;  getK(param);
    SparseMatrix<double, 0, int> getK_by_eps(int pos, double eps) {
        VectorXd tmp = theta_K;
        tmp(pos) = tmp(pos) + eps;
        return getK(tmp);
    }

    SparseMatrix<double, 0, int> get_dK_by_eps(int index, int pos, double eps) {
        VectorXd tmp = theta_K;
        tmp(pos) = tmp(pos) + eps;
        return get_dK(index, tmp);
    }

    /* 4 for optimizer */
    const VectorXd get_parameter() const;
    const VectorXd get_grad();
    void           set_parameter(const VectorXd&);
    void           finishOpt(int i) {fix_flag[i] = 0; }


    // deprecated
    // virtual double function_kappa(double eps);
    virtual double function_K(VectorXd& parameter);

    // used for general case
    virtual double function_K(SparseMatrix<double>& K);
    virtual VectorXd numerical_grad(); // given eps

    // update the trace value
    void compute_trace() {
        if (W_size != V_size) return;

// compute trace
// auto timer_trace = std::chrono::steady_clock::now();

        SparseMatrix<double> K = getK(theta_K);
        SparseMatrix<double> dK = get_dK_by_index(0);
        if (!numer_grad) {
            if (!symmetricK) {
                lu_solver_K.computeKTK(K);
                trace = lu_solver_K.trace(dK);
            } else {
                chol_solver_K.compute(K);
                trace = chol_solver_K.trace(dK);
            }
        }
// std::cout << "eps K  in 1= " << K << std::endl;
// std::cout << "eps dK in 1= " << dK << std::endl;
// std::cout << "trace  in 1= " << trace << std::endl;
        // else {
        //     if (!symmetricK) {
        //         // BiCG solver
        //         throw("Not implemented yet");
        //     } else {
        //         SparseMatrix<double, RowMajor> KK = K;
        //         CG_solver_K.compute(K);
        //         trace = CG_solver_K.trace(M);
        //     }
        // }

// std::cout << "trace ====== " << trace << std::endl;
// std::cout << "time for the trace (ms): " << since(timer_trace).count() << std::endl;

        // update trace_eps if using hessian
        if ((!numer_grad) && (use_precond)) {
            SparseMatrix<double> K = getK_by_eps(0, eps);
            SparseMatrix<double> dK = get_dK_by_eps(0, 0, eps);
            SparseMatrix<double> M = dK;

            if (!symmetricK) {
                lu_solver_K.computeKTK(K);
                trace_eps = lu_solver_K.trace(M);
            } else {
                chol_solver_K.compute(K);
                trace_eps = chol_solver_K.trace(M);
            }

            // else {
            //     if (!symmetricK) {
            //         // BiCG solver
            //         throw("Not implemented yet");
            //     } else {
            //         CG_solver_K.compute(K);
            //         trace_eps = CG_solver_K.trace(M);
            //     }
            // }
        }
    };

    virtual VectorXd grad_theta_K() { return numerical_grad(); }
    virtual VectorXd grad_theta_mu();
    virtual VectorXd grad_theta_sigma();
    virtual VectorXd grad_theta_sigma_normal(); // grad of sig. only for normal noise

    // mean of grad for each var
    double grad_theta_nu() {
        double grad = 0;
        for (int i=0; i < n_rep; i++)
            grad += vars[i].grad_log_nu();
        return grad / n_rep;
    }

    // Output
    // virtual Rcpp::List get_estimates() const=0;
    Rcpp::List output() const;

    void record_traj() {
        for (int i=0; i < theta_K.size(); i++)
            theta_K_traj[i].push_back(theta_K(i));
        for (int i=0; i < theta_mu.size(); i++)
            theta_mu_traj[i].push_back(theta_mu(i));
        for (int i=0; i < theta_sigma.size(); i++)
            theta_sigma_traj[i].push_back(theta_sigma(i));
        if (noise_type=="normal_nig")
            for (int i=0; i < theta_sigma_normal.size(); i++)
                theta_sigma_normal_traj[i].push_back(theta_sigma_normal(i));
        nu_traj.push_back(vars[0].get_nu());
    }

    // stop estimate theta_K if converged
    // void check_converge(vector<bool>& converge) {
    //     vector<bool> K_converge (converge.begin(), converge.begin() + n_theta_K);
    //     // if all K_converge fix
    //     bool all_converge = true;
    //     for (int i=0; i<n_theta_K; i++)
    //         if (!K_converge[i]) all_converge = false;
    //     if (all_converge) fix_flag[latent_fix_theta_K] = 1;
    // }
};

/*    Optimizer related    */
inline const VectorXd Latent::get_parameter() const {
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

if (debug) std::cout << "parameter= " << parameter << std::endl;
// if (debug) std::cout << "End latent get parameter"<< std::endl;
    return parameter;
}

inline const VectorXd Latent::get_grad() {
if (debug) std::cout << "Start latent gradient"<< std::endl;
// auto grad1 = std::chrono::steady_clock::now();
    VectorXd grad = VectorXd::Zero(n_params);

    if (noise_type == "normal") {
        if (!fix_flag[latent_fix_theta_K])
            grad.segment(0, n_theta_K) = grad_theta_K();
if (debug) std::cout << "after K"<< std::endl;
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
    std::cout << "gradient= " << grad << std::endl;
    // std::cout << "one latent gradient time " << since(grad1).count() << std::endl;
}
if (debug) std::cout << "finish latent gradient"<< std::endl;
    return grad;
}

inline void Latent::set_parameter(const VectorXd& theta) {
// if (debug) std::cout << "Start latent set parameter"<< std::endl;
    if (noise_type == "normal") {
        theta_K  = theta.segment(0, n_theta_K);
        theta_sigma = theta.segment(n_theta_K, n_theta_sigma);
        sigma = (B_sigma * theta_sigma).array().exp();
    } else {
        // nig, gal and normal+nig
        theta_K  = theta.segment(0, n_theta_K);
        theta_mu = theta.segment(n_theta_K, n_theta_mu);
        theta_sigma = theta.segment(n_theta_K+n_theta_mu, n_theta_sigma);
        double log_nu = (theta(n_theta_K+n_theta_mu+n_theta_sigma));
        for (int i=0; i < n_rep; i++) vars[i].set_log_nu(log_nu); // for each replicate

        // update
        mu = (B_mu * theta_mu);
        sigma = (B_sigma * theta_sigma).array().exp();
        if (noise_type == "normal_nig") {
            theta_sigma_normal = theta.segment(n_theta_K+n_theta_mu+n_theta_sigma+1, n_theta_sigma_normal);
            sigma_normal = (B_sigma_normal * theta_sigma_normal).array().exp();
        }
    }

    // update K, K_rep, ...
    update_each_iter();

    // record
    record_traj();
}

// subclasses
class AR : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
    bool is_rw; // if this is actually a rw model
public:
    AR(const Rcpp::List& model_list, unsigned long seed, bool is_rw);

    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
    VectorXd grad_theta_K();
    void update_each_iter();
    void update_num_dK();
    double th2a(double th) const {return (-1 + 2*exp(th) / (1+exp(th)));}
    double a2th(double k) const {return (log((-1-k)/(-1+k)));}

    VectorXd unbound_to_bound_K(const VectorXd& theta_K) const {
        VectorXd alpha (1);
        alpha(0) = th2a(theta_K(0));
        return alpha;
    }
    VectorXd bound_to_unbound_K(const VectorXd& alpha) const {
        VectorXd theta_K (1);
        theta_K(0) = a2th(alpha(0));
        return theta_K;
    }
    // Rcpp::List get_estimates() const {
    //     return Rcpp::List::create(
    //         Rcpp::Named("alpha")        = theta_K(0),
    //         Rcpp::Named("theta.mu")     = theta_mu,
    //         Rcpp::Named("theta.sigma")  = theta_sigma,
    //         Rcpp::Named("theta.noise")  = var.get_nu()
    //     );
    // }
};


class Matern : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    VectorXd Cdiag;
public:
    Matern(const Rcpp::List& model_list, unsigned long seed);
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
    // Rcpp::List get_estimates() const {
    //     return Rcpp::List::create(
    //         Rcpp::Named("kappa")        = theta_K(0),
    //         Rcpp::Named("theta.mu")     = theta_mu,
    //         Rcpp::Named("theta.sigma")  = theta_sigma,
    //         Rcpp::Named("theta.noise")  = var.get_nu()
    //     );
    // }
};

class Matern_ns : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    MatrixXd Bkappa;
    VectorXd Cdiag;
public:
    Matern_ns(const Rcpp::List& model_list, unsigned long seed);
    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
    VectorXd grad_theta_K();
    void update_each_iter();

    // Rcpp::List get_estimates() const {
    //     return Rcpp::List::create(
    //         Rcpp::Named("theta.kappa") = theta_K,
    //         Rcpp::Named("theta.mu")    = theta_mu,
    //         Rcpp::Named("theta.sigma") = theta_sigma,
    //         Rcpp::Named("theta.noise") = var.get_nu()
    //     );
    // }
};

// for initialize Latent models
class LatentFactory {
public:
  static std::unique_ptr<Latent> create(const std::string& model_type, const Rcpp::List& latent_in, int latent_seed) {
    int n_theta_K = Rcpp::as<int> (latent_in["n_theta_K"]);

    if (model_type == "ar1") {
      return std::make_unique<AR>(latent_in, latent_seed, false);
    } else if (model_type == "rw1") {
      return std::make_unique<AR>(latent_in, latent_seed, true);
    } else if (model_type == "matern" && n_theta_K > 1) {
      return std::make_unique<Matern_ns>(latent_in, latent_seed);
    } else if (model_type == "matern" && n_theta_K == 1) {
      return std::make_unique<Matern>(latent_in, latent_seed);
    } else {
      throw std::runtime_error("Unknown model.");
    }
  }
};

#endif