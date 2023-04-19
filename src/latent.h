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
#include <memory>

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
public:
    bool symmetricK {false};
protected:
    std::mt19937 latent_rng;
    string model_type, noise_type;
    bool debug;

    int n_rep, W_size, V_size, n_params, n_var {1}; // n_params=n_theta_K + n_theta_mu + n_theta_sigma + n_var

    // operator K related
    VectorXd theta_K;
    int n_theta_K;

    bool use_num_dK {false};
    SparseMatrix<double, 0, int> K; // size = V_size * W_size
    vector<SparseMatrix<double, 0, int>> dK;
    vector<double> trace;
    bool zero_trace;

    SparseMatrix<double, 0, int> K_rep;  // n_rep of size of K (for sampling W)

    bool fix_flag[LATENT_FIX_FLAG_SIZE] {0};

    bool use_precond {false}, numer_grad {false};

    // mu and sigma, and sigma_normal (special case when using nig_normal case)
    MatrixXd B_mu, B_sigma, B_sigma_normal;
    VectorXd theta_mu, theta_sigma, theta_sigma_normal;

    // mu = Bmu * theta_mu
    // sigma = exp(Bsigma * theta_sigma)
    VectorXd mu, sigma, sigma_normal;
    int n_theta_mu, n_theta_sigma, n_theta_sigma_normal;

    // for numerical gradient.
    double eps;

    // VectorXd W, prevW, h;
    VectorXd h;
    SparseMatrix<double,0,int> A;

    vector<VectorXd> Ws, prevWs;
    vector<Var> vars;

    // solver
    cholesky_solver chol_solver_K;
    lu_sparse_solver lu_solver_K;
    bool use_iter_solver {false};
    // iterative_solver CG_solver_K;

    cholesky_solver solver_Q; // Q = KT diag(1/SV) K
public:
    Latent(const Rcpp::List&, unsigned long seed);
    virtual ~Latent() {}

    virtual void update_each_iter();

    /*  1 Model itself   */
    int get_W_size() const                  {return n_rep * W_size; }
    int get_V_size() const                  {return n_rep * V_size; }
    int get_n_params() const                {return n_params; }
    int get_n_theta_K() const               {return n_theta_K; }
    void set_theta_K(const VectorXd& theta) {theta_K = theta; }

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
            // V-h or V-1
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

    virtual void sample_cond_V() {
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
    virtual const VectorXd get_grad();
    void           set_parameter(const VectorXd&);
    void           finishOpt(int i) {fix_flag[i] = 0; }

    // used for general case
    virtual double function_K(VectorXd& parameter);
    virtual double function_K(SparseMatrix<double>& K);

    virtual VectorXd grad_theta_K();
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

    Rcpp::List output() const;
};


// subclasses
enum Type {ar, rw, ou, matern_ns};
class AR : public Latent {
private:
    SparseMatrix<double, 0, int> G, C;
    Type type;
    MatrixXd B_K;
public:
    AR(const Rcpp::List& model_list, unsigned long seed, Type type);
    SparseMatrix<double> getK(const VectorXd& alpha) const override;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const override;

    double th2a(double th) const {return (-1 + 2*exp(th) / (1+exp(th)));}
    double a2th(double k) const {return (log((-1-k)/(-1+k)));}
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
};

class Matern_ns : public Latent {
private:
    Type type;
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    MatrixXd Bkappa;
    VectorXd Cdiag;
public:
    Matern_ns(const Rcpp::List& model_list, unsigned long seed, Type type);
    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
};

class Tensor_prod : public Latent {
private:
  int n_theta_l, n_theta_r, V_size_l, V_size_r, W_size_l, W_size_r;
  std::unique_ptr<Latent> left, right;
public:
  Tensor_prod(const Rcpp::List& model_list, unsigned long seed);
  SparseMatrix<double> getK(const VectorXd& alpha) const;
  SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
};

class Iid : public Latent {
private:
    SparseMatrix<double, 0, int> I;
public:
  Iid(const Rcpp::List& model_list, unsigned long seed);
  SparseMatrix<double> getK(const VectorXd& alpha) const;
  SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
};

// ---- Structure for random effects ----
// U|V ~ N(0, Sigma)
class Randeff : public Latent{
  private:
    std::mt19937 randeff_rng;
    int n_repl;
  public:
    Randeff(const Rcpp::List& R_randeff, unsigned long seed);
    SparseMatrix<double> getK(const VectorXd& theta_K) const;
    MatrixXd get_dK_dense(int index, const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
    VectorXd grad_theta_K();
    void sample_cond_V() override;
    void update_each_iter();
};

// Bivar
class Bivar : public Latent {
private:
    std::unique_ptr<Latent> m1, m2;
    int n; // dim of K1 and K2 (same)
    bool share_param;
public:
    Bivar(const Rcpp::List& model_list, unsigned long seed);
    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;

    Matrix2d getD(double, double) const;
    Matrix2d get_dD_theta(double, double) const;
    Matrix2d get_dD_rho(double, double) const;
    Matrix2d get_dD2_theta(double, double) const;
    Matrix2d get_dD2_rho(double, double) const;
};

// for initialize Latent models
class LatentFactory {
public:
  static std::unique_ptr<Latent> create(const Rcpp::List& latent_in, int latent_seed) {
    int n_theta_K = Rcpp::as<int> (latent_in["n_theta_K"]);
    string model_type = Rcpp::as<string> (latent_in["model"]);

    if (model_type == "tp") {
      return std::make_unique<Tensor_prod>(latent_in, latent_seed);
    } else if (model_type == "ar1") {
      return std::make_unique<AR>(latent_in, latent_seed, Type::ar);
    } else if (model_type == "rw") {
      return std::make_unique<AR>(latent_in, latent_seed, Type::rw);
    } else if (model_type == "ou") {
      return std::make_unique<Matern_ns>(latent_in, latent_seed, Type::ou);
    } else if (model_type == "matern" && n_theta_K > 1) {
      return std::make_unique<Matern_ns>(latent_in, latent_seed, Type::matern_ns);
    } else if (model_type == "matern" && n_theta_K == 1) {
      return std::make_unique<Matern>(latent_in, latent_seed);
    } else if (model_type == "iid") {
      return std::make_unique<Iid>(latent_in, latent_seed);
    } else if (model_type == "re") {
      return std::make_unique<Randeff>(latent_in, latent_seed);
    } else if (model_type == "bv") {
      return std::make_unique<Bivar>(latent_in, latent_seed);
    } else {
      throw std::runtime_error("Unknown model.");
    }
  }
};

#endif