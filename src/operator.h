#ifndef NGME_OPERATOR_H
#define NGME_OPERATOR_H

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
enum Type {ar, rw, ou, matern_ns};

class Operator {
protected:
    VectorXd h;
    int n_theta_K;
    bool zero_trace, symmetric;
public:
    Operator(const Rcpp::List& operator_list) :
        h (Rcpp::as<VectorXd> (operator_list["h"])),
        n_theta_K (Rcpp::as<int> (operator_list["n_theta_K"])),
        zero_trace (Rcpp::as<bool> (operator_list["zero_trace"])),
        symmetric (Rcpp::as<bool> (operator_list["symmetric"]))
    {}

    int get_n_theta_K() const {return n_theta_K;}
    const VectorXd& get_h() const {return h;}
    bool is_symmetric() const {return symmetric;}

    virtual SparseMatrix<double, 0, int> getK(const VectorXd&) const = 0;
    virtual SparseMatrix<double, 0, int> get_dK(int index, const VectorXd&) const = 0;
};

class AR : public Operator {
private:
    SparseMatrix<double, 0, int> G, C;
    Type type;
    MatrixXd B_K;
public:
    AR(const Rcpp::List&, Type type);

    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;

    double th2a(double th) const {return (-1 + 2*exp(th) / (1+exp(th)));}
    double a2th(double k) const {return (log((-1-k)/(-1+k)));}
};

class Matern : public Operator {
private:
    SparseMatrix<double, 0, int> G, C;
    int alpha;
    VectorXd Cdiag;
public:
    Matern(const Rcpp::List&);

    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
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

    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
};

class Tensor_prod : public Operator {
private:
  int n_theta_l, n_theta_r, V_size_l, V_size_r, W_size_l, W_size_r;
  std::unique_ptr<Operator> left, right;
public:
  Tensor_prod(const Rcpp::List&);

  SparseMatrix<double> getK(const VectorXd& alpha) const;
  SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
};

class Iid : public Operator {
private:
    SparseMatrix<double, 0, int> I;
public:
  Iid(const Rcpp::List&);

  SparseMatrix<double> getK(const VectorXd& alpha) const;
  SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
};

// ---- Structure for random effects ----
// U|V ~ N(0, Sigma)
class Randeff : public Operator{
  private:
    std::mt19937 randeff_rng;
    int n_repl;
    int W_size;
  public:
    Randeff(const Rcpp::List&);

    SparseMatrix<double> getK(const VectorXd& theta_K) const;
    MatrixXd get_dK_dense(int index, const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;
    // VectorXd grad_theta_K();
    // void sample_cond_V() override;
    // void update_each_iter();
};

// Bivar
class Bivar : public Operator {
private:
    std::unique_ptr<Operator> m1, m2;
    int n; // dim of K1 and K2 (same)
    bool share_param;
public:
    Bivar(const Rcpp::List&);

    SparseMatrix<double> getK(const VectorXd& alpha) const;
    SparseMatrix<double> get_dK(int index, const VectorXd& alpha) const;

    Matrix2d getD(double, double) const;
    Matrix2d get_dD_theta(double, double) const;
    Matrix2d get_dD_rho(double, double) const;
    Matrix2d get_dD2_theta(double, double) const;
    Matrix2d get_dD2_rho(double, double) const;
};

// for initialize Latent models
class OperatorFactory {
public:
  static std::unique_ptr<Operator> create(
    const Rcpp::List& operator_in
  ) {
    string model_type = Rcpp::as<string> (operator_in["model"]);
    VectorXd theta_K = Rcpp::as<VectorXd> (operator_in["theta_K"]);
    int n_theta_K = theta_K.size();

    if (model_type == "tp") {
      return std::make_unique<Tensor_prod>(operator_in);
    } else if (model_type == "ar1") {
      return std::make_unique<AR>(operator_in, Type::ar);
    } else if (model_type == "rw") {
      return std::make_unique<AR>(operator_in, Type::rw);
    } else if (model_type == "ou") {
      return std::make_unique<Matern_ns>(operator_in, Type::ou);
    } else if (model_type == "matern" && n_theta_K > 1) {
      return std::make_unique<Matern_ns>(operator_in, Type::matern_ns);
    } else if (model_type == "matern" && n_theta_K == 1) {
      return std::make_unique<Matern>(operator_in);
    } else if (model_type == "iid") {
      return std::make_unique<Iid>(operator_in);
    } else if (model_type == "re") {
      return std::make_unique<Randeff>(operator_in);
    } else if (model_type == "bv") {
      return std::make_unique<Bivar>(operator_in);
    } else {
      throw std::runtime_error("Unknown model.");
    }
  }
};

#endif