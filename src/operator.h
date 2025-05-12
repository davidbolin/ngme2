#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

// #include<Eigen/IterativeLinearSolvers>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <memory>

#include "include/timer.h"
#include "include/solver.h"

using std::exp;
using std::log;
using std::pow;
using std::string;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::VectorXd;
using std::vector;


// subclasses
enum Type {ar, ou, matern_ns};

class Operator {
protected:
    VectorXd h;
    int n_theta_K;
    bool zero_trace, symmetric, generic;

    SparseMatrix<double> K;
    vector<SparseMatrix<double>> dK;
public:
    Operator(const Rcpp::List& operator_list) :
        h (Rcpp::as<VectorXd> (operator_list["h"])),
        n_theta_K (Rcpp::as<int> (operator_list["n_theta_K"])),
        zero_trace (Rcpp::as<bool> (operator_list["zero_trace"])),
        symmetric (Rcpp::as<bool> (operator_list["symmetric"])),
        generic (Rcpp::as<bool> (operator_list["generic"])),
        K (Rcpp::as<SparseMatrix<double>> (operator_list["K"])),
        dK (n_theta_K)
    {
      // initial dK
      for (int i=0; i<n_theta_K; i++) {
        dK[i].resize(h.size(), h.size());
        dK[i].setZero();
      }
    }
    virtual ~Operator() = default;

    int get_n_theta_K() const {return n_theta_K;}
    const VectorXd& get_h() const {return h;}
    bool is_symmetric() const {return symmetric;}
    bool is_zero_trace() const {return zero_trace;}

    const SparseMatrix<double>& getK() const {return K;}
    const vector<SparseMatrix<double>>& get_dK() const {return dK;}

    virtual void update_K(const VectorXd&) = 0;
    virtual void update_dK(const VectorXd&) = 0;
};

class Matern : public Operator {
private:
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    VectorXd Cdiag;
public:
    Matern(const Rcpp::List&);

    void update_K(const VectorXd&);
    void update_dK(const VectorXd&);

    int get_alpha() const {return alpha;}
};

class Matern_ns : public Operator {
private:
    Type type;
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    MatrixXd Bkappa;
    VectorXd Cdiag;
public:
    Matern_ns(const Rcpp::List&, Type);

    void update_K(const VectorXd&);
    void update_dK(const VectorXd&);
};

class Tensor_prod : public Operator {
private:
  std::shared_ptr<Operator> first, second;
  int n_theta_1, n_theta_2;
public:
  Tensor_prod(const Rcpp::List&);

  void update_K(const VectorXd&);
  void update_dK(const VectorXd&);
};

class Spacetime : public Operator {
private:
  VectorXd Ct_diag, Cs_diag; // not used
  SparseMatrix<double, 0, int> BtCs, Gs, Ct, Cs, Bx, By, S, Bs, Hxx, Hyy, Hxy, Hyx;
  // MatrixXd B_gamma_x, B_gamma_y;
  Rcpp::List B_gamma_x_list_input, B_gamma_y_list_input;
  VectorXd theta_gamma_x, theta_gamma_y;
  int n_theta_gamma_x, n_theta_gamma_y;
  double lambda, alpha;
  string method; // galerkin, backward Euler
  bool stabilization, fix_gamma;
  int nt;
  std::vector<MatrixXd> B_gamma_x_list, B_gamma_y_list;
public:
  Spacetime(const Rcpp::List&);

  void update_K(const VectorXd&);
  void update_dK(const VectorXd&);
};

// Bivar
class Bivar : public Operator {
private:
    std::shared_ptr<Operator> first, second;
    int n_theta_1, n_theta_2;
    int n; // dim of K1 and K2 (same)
    bool share_param, fix_bv_theta;
public:
    Bivar(const Rcpp::List&);

    void update_K(const VectorXd&);
    void update_dK(const VectorXd&);

    Matrix2d getD(double, double) const;
    Matrix2d get_dD_theta(double, double) const;
    Matrix2d get_dD_rho(double, double) const;
    Matrix2d get_dD2_theta(double, double) const;
    Matrix2d get_dD2_rho(double, double) const;
};

class Generic : public Operator {
private:
    vector<SparseMatrix<double, 0, int>> matrices;
    MatrixXd idx_mat;
    vector<string> trans;
public:
    Generic(const Rcpp::List&);

    void update_K(const VectorXd& theta_K);
    void update_dK(const VectorXd& theta_K);

    // VectorXd param_trans_fun(const VectorXd& theta_K, const string& trans_type) const;
    VectorXd compute_coef(const VectorXd& theta_K, const MatrixXd& idx_mat, const vector<string>& trans) const;
};

// Bivar_normal
// class Bivar_normal : public Operator {
// private:
//     std::shared_ptr<Operator> first, second;
//     int n_theta_1, n_theta_2;
//     int n; // dim of K1 and K2 (same)
//     bool share_param, fix_bv_theta;
// public:
//     Bivar_normal(const Rcpp::List&);

//     void update_K(const VectorXd&);
//     void update_dK(const VectorXd&);

//     Matrix2d getD(double, double) const;
//     Matrix2d get_dD_theta(double, double) const;
//     Matrix2d get_dD_rho(double, double) const;
//     Matrix2d get_dD2_theta(double, double) const;
//     Matrix2d get_dD2_rho(double, double) const;
// };


// Bivar_normal_ope (theta=0)
class Bivar_normal_ope : public Operator {
private:
    std::shared_ptr<Operator> first, second;
    int n_theta_1, n_theta_2;
    int n; // dim of K1 and K2 (same)
    bool share_param, fix_bv_theta;
public:
    Bivar_normal_ope(const Rcpp::List&);

    void update_K(const VectorXd&);
    void update_dK(const VectorXd&);

    Matrix2d getD(double, double) const;
    Matrix2d get_dD_theta(double, double) const;
    Matrix2d get_dD_rho(double, double) const;
    Matrix2d get_dD2_theta(double, double) const;
    Matrix2d get_dD2_rho(double, double) const;
};


class bv_matern_normal : public Operator {
private:
    std::shared_ptr<Matern> first, second;
    int n_theta_1, n_theta_2;
    int n; // dim of K1 and K2 (same)
    bool share_param, fix_bv_theta;
    double dim, alpha1, alpha2, nu1, nu2;
public:
    bv_matern_normal(const Rcpp::List&);

    void update_K(const VectorXd&);
    void update_dK(const VectorXd&);

    Matrix2d getD(double, double) const;
    Matrix2d get_dD_theta(double, double) const;
    Matrix2d get_dD_rho(double, double) const;
    Matrix2d get_dD2_theta(double, double) const;
    Matrix2d get_dD2_rho(double, double) const;
};

class bv_matern_nig : public Operator {
private:
    std::shared_ptr<Matern> first, second;
    int n_theta_1, n_theta_2;
    int n; // dim of K1 and K2 (same)
    bool share_param, fix_bv_theta;
    double dim, alpha1, alpha2, nu1, nu2, bv_theta;
public:
    bv_matern_nig(const Rcpp::List&);

    void update_K(const VectorXd&);
    void update_dK(const VectorXd&);

    Matrix2d getD(double, double) const;
    Matrix2d get_dD_theta(double, double) const;
    Matrix2d get_dD_rho(double, double) const;
    Matrix2d get_dD2_theta(double, double) const;
    Matrix2d get_dD2_rho(double, double) const;
};


// ---- Structure for random effects ----
// U|V ~ N(0, Sigma)
class Randeff : public Operator{
  private:
    int n_reff;
  public:
    Randeff(const Rcpp::List&);

    void update_K(const VectorXd& theta_K);
    void update_dK(const VectorXd& theta_K);
};

// for initialize Latent models
class OperatorFactory {
public:
  static std::shared_ptr<Operator> create(
    const Rcpp::List& operator_in
  ) {
    string model_type = Rcpp::as<string> (operator_in["model"]);
    VectorXd theta_K = Rcpp::as<VectorXd> (operator_in["theta_K"]);
    bool generic = Rcpp::as<bool> (operator_in["generic"]);
    int n_theta_K = theta_K.size();

    if (model_type == "generic") {
      return std::make_shared<Generic>(operator_in);
    } else if (generic) {
      // using the generic structure to build the operator
      return std::make_shared<Generic>(operator_in["generic_operator"]);
    } else if (model_type == "tp") {
      return std::make_shared<Tensor_prod>(operator_in);
    } else if (model_type == "spacetime") {
      return std::make_shared<Spacetime>(operator_in);
    } else if (model_type == "ou") {
      return std::make_shared<Matern_ns>(operator_in, Type::ou);
    } else if (model_type == "matern" && n_theta_K > 1) {
      return std::make_shared<Matern_ns>(operator_in, Type::matern_ns);
    } else if (model_type == "matern" && n_theta_K == 1) {
      return std::make_shared<Matern>(operator_in);
    } else if (model_type == "re") {
      return std::make_shared<Randeff>(operator_in);
    } else if (model_type == "bv") {
      return std::make_shared<Bivar>(operator_in);
    } else if (model_type == "bv_normal") {
      return std::make_shared<Bivar_normal_ope>(operator_in);
    } else if (model_type == "bv_matern_normal") {
      return std::make_shared<bv_matern_normal>(operator_in);
    } else if (model_type == "bv_matern_nig") {
      return std::make_shared<bv_matern_nig>(operator_in);
    } else {
      throw std::runtime_error("Unknown model.");
    }
  }
};

#endif
