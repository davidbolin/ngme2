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
enum Type {ar, ou, matern_ns};

class Operator {
protected:
    VectorXd h;
    int n_theta_K;
    bool zero_trace, symmetric;

    SparseMatrix<double> K;
    vector<SparseMatrix<double>> dK;
public:
    Operator(const Rcpp::List& operator_list) :
        h (Rcpp::as<VectorXd> (operator_list["h"])),
        n_theta_K (Rcpp::as<int> (operator_list["n_theta_K"])),
        zero_trace (Rcpp::as<bool> (operator_list["zero_trace"])),
        symmetric (Rcpp::as<bool> (operator_list["symmetric"])),
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

class AR : public Operator {
private:
    SparseMatrix<double, 0, int> G, C;
    Type type;
    MatrixXd B_K;
public:
    AR(const Rcpp::List&);

    void update_K(const VectorXd& alpha);
    void update_dK(const VectorXd& alpha);

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
  SparseMatrix<double, 0, int> BtCs,   Gs, Bs, Ct, Cs, S;
  double lambda, alpha;
  string method; // galerkin, backward Euler
  bool stabilization;
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




// notice dK is of size 0
class Iid : public Operator {
public:
  Iid(const Rcpp::List& operator_list):
    Operator(operator_list)
  {}

  void update_K(const VectorXd& alpha) {};
  void update_dK(const VectorXd& alpha) {};
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
    int n_theta_K = theta_K.size();

    if (model_type == "tp") {
      return std::make_shared<Tensor_prod>(operator_in);
    } else if (model_type == "ar1") {
      return std::make_shared<AR>(operator_in);
    } else if (model_type == "spacetime") {
      return std::make_shared<Spacetime>(operator_in);
    } else if (model_type == "ou") {
      return std::make_shared<Matern_ns>(operator_in, Type::ou);
    } else if (model_type == "matern" && n_theta_K > 1) {
      return std::make_shared<Matern_ns>(operator_in, Type::matern_ns);
    } else if (model_type == "matern" && n_theta_K == 1) {
      return std::make_shared<Matern>(operator_in);
    } else if (model_type == "iid" || model_type == "rw1" || model_type == "rw2") {
      return std::make_shared<Iid>(operator_in);
    } else if (model_type == "re") {
      return std::make_shared<Randeff>(operator_in);
    } else if (model_type == "bv") {
      return std::make_shared<Bivar>(operator_in);    
    } else if (model_type == "bv_matern_normal") {
      return std::make_shared<bv_matern_normal>(operator_in);
    } else if (model_type == "bv_normal") {
      return std::make_shared<Bivar_normal_ope>(operator_in);
    } else {
      throw std::runtime_error("Unknown model.");
    }
  }
};

#endif