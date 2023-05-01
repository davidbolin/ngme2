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
#include "operator.h"

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
    std::unique_ptr<Operator> ope;

    VectorXd theta_K;
    int n_theta_K;
    bool symmetricK;

    int n_rep, W_size, V_size, n_params, n_var {1}; // n_params=n_theta_K + n_theta_mu + n_theta_sigma + n_var

    bool use_num_dK {false};
    SparseMatrix<double, 0, int> K; // size = V_size * W_size
    vector<SparseMatrix<double, 0, int>> dK;
    vector<double> trace;
    bool zero_trace;

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
    cholesky_solver solver_Q; // Q = KT diag(1/SV) K
public:
    Latent(const Rcpp::List&, unsigned long seed);
    ~Latent() {}

    void update_each_iter();

    /*  1 Model itself   */
    int get_W_size() const             {return n_rep * W_size; }
    int get_V_size() const             {return n_rep * V_size; }
    int get_n_params() const           {return n_params; }
    int get_n_theta_K() const          {return n_theta_K; }
    const VectorXd& get_theta_K() const {return theta_K; }

    SparseMatrix<double, 0, int>& getA()  {return A;}

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
        return K;
    }

    // param(pos) += eps;  getK(param);
    // SparseMatrix<double, 0, int> getK_by_eps(int pos, double eps) {
    //     VectorXd tmp = get_theta_K();
    //     tmp(pos) = tmp(pos) + eps;
    //     return ope->getK(tmp);
    // }

    // SparseMatrix<double, 0, int> get_dK_by_eps(int index, int pos, double eps) {
    //     VectorXd tmp = get_theta_K();
    //     tmp(pos) = tmp(pos) + eps;
    //     return get_dK(index, tmp);
    // }

    /* 4 for optimizer */
    const VectorXd get_parameter() const;
    const VectorXd get_grad();
    void           set_parameter(const VectorXd&);
    void           finishOpt(int i) {fix_flag[i] = 0; }

    // used for general case
    double function_K(VectorXd& parameter);
    double function_K(SparseMatrix<double>& K);

    VectorXd grad_theta_K();
    VectorXd grad_theta_mu();
    VectorXd grad_theta_sigma();
    VectorXd grad_theta_sigma_normal(); // grad of sig. only for normal noise

    // mean of grad for each var
    double grad_theta_nu() {
        double grad = 0;
        for (int i=0; i < n_rep; i++)
            grad += vars[i].grad_log_nu();
        return grad / n_rep;
    }

    Rcpp::List output() const;
};

#endif