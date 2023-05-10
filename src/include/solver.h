#ifndef __Solver__Solver__
#define __Solver__Solver__
#include <iostream>
#include <string.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include "MatrixAlgebra.h"
#include <Rcpp.h>

class solver
{
private:
public:
  int n;
  solver(const solver &){};
  solver(){};
  virtual ~solver(){};
  virtual void init(int, int, int, double) = 0;
  virtual void initFromList(int, Rcpp::List const &) = 0;
  virtual inline void analyze(const Eigen::SparseMatrix<double, 0, int> &) = 0;
  virtual void compute(const Eigen::SparseMatrix<double, 0, int> &) = 0;
  virtual void compute(const Eigen::MatrixXd &) { std::cout << "compute not implimented for MatrixXd\n"; };
  virtual inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &) = 0;
  virtual double trace(const Eigen::MatrixXd &) = 0;
  virtual double trace(const Eigen::SparseMatrix<double, 0, int> &) = 0;
  virtual double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &) = 0;
  virtual inline double logdet() { return 0.0; }
  virtual inline Eigen::VectorXd Qinv_diag()
  {
    Eigen::VectorXd a;
    a.setZero(1);
    return a;
  }
  virtual Eigen::VectorXd rMVN(Eigen::VectorXd &, Eigen::VectorXd &) = 0;
  virtual SparseMatrix<double, 0, int> return_Qinv()
  {
    SparseMatrix<double, 0, int> a;
    std::cout << "return_Qinv not implimented for MatrixXd\n";
    return a;
  }
};

// class iterative_solver : public virtual solver
// {
// private:
//   int N;
//   Eigen::MatrixXd U, QU, MQU, MU;
//   bool QU_computed;
//   Eigen::ConjugateGradient<Eigen::SparseMatrix<double, 0, int>, Lower|Upper > R;

// public:
//   ~iterative_solver(){};
//   void init(int, int, int, double);
//   void initFromList(int, Rcpp::List const &);
//   inline void analyze(Eigen::SparseMatrix<double, 0, int> &M) { R.analyzePattern(M); }
//   void compute(Eigen::SparseMatrix<double, 0, int> &);
//   inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return R.solveWithGuess(v, x); }
//   double trace(Eigen::MatrixXd &);
//   double trace(Eigen::SparseMatrix<double, 0, int> &);
//   double trace2(SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);
//   Eigen::VectorXd rMVN(Eigen::VectorXd &mu, Eigen::VectorXd &z)
//   {
//     std::cout << "rMVN not implimented\n";
//     Eigen::VectorXd X;
//     return X;
//   }
// };


class cholesky_solver : public virtual solver
{
private:
  bool Qi_computed;
  double ld;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int> > R;
  Eigen::SparseMatrix<double, 0, int> Qi;
  void set_ld();

public:
  cholesky_solver(const cholesky_solver &){};
  cholesky_solver(){};
  ~cholesky_solver(){};
  void init(int, int, int, double);
  void initFromList(int, Rcpp::List const &);
  inline void analyze(const Eigen::SparseMatrix<double, 0, int> &M) { R.analyzePattern(M); }
  void compute(const Eigen::SparseMatrix<double, 0, int> &);
  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return R.solve(v); }
  inline Eigen::VectorXd solve(const Eigen::VectorXd &v)               { return R.solve(v); }
  double trace(const Eigen::MatrixXd &);
  double trace(const Eigen::SparseMatrix<double, 0, int> &);
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);
  inline double logdet()
  {
    set_ld();
    return ld;
  }
  VectorXd Qinv_diag();
  Eigen::VectorXd rMVN(Eigen::VectorXd &, Eigen::VectorXd &);
  SparseMatrix<double, 0, int> return_Qinv();
};



class lu_solver : public virtual solver
{
private:
  Eigen::MatrixXd Kinv;
  Eigen::FullPivLU<Eigen::MatrixXd> LU_K;
  int Kinv_computed;

public:
  void init(int, int, int, double);
  void initFromList(int, Rcpp::List const &);
  void analyze(const Eigen::SparseMatrix<double, 0, int> &);
  void compute(const Eigen::SparseMatrix<double, 0, int> &);
  void compute(const Eigen::MatrixXd &);
  Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &);
  double trace(const Eigen::MatrixXd &);
  double trace(const Eigen::SparseMatrix<double, 0, int> &);
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);
  double logdet();
  Eigen::VectorXd rMVN(Eigen::VectorXd &, Eigen::VectorXd &)
  {
    std::cout << "lu_solver:rMVN not implimented\n";
    throw;
  };
};

class lu_sparse_solver : public virtual solver
{
private:
  int n;
  Eigen::SparseMatrix<double, 0, int> KKtinv;
  Eigen::SparseMatrix<double, 0, int> K;
  Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int> > LU_K;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int> > L_KKt;
  int KKtinv_computed;

public:
  void init(int, int, int, double);
  void initFromList(int, Rcpp::List const &);
  void analyze(const Eigen::SparseMatrix<double, 0, int> &);
  void compute(const Eigen::SparseMatrix<double, 0, int> &);
  void computeKTK(Eigen::SparseMatrix<double, 0, int> &);
  double trace(const Eigen::MatrixXd &);
  double trace(const Eigen::SparseMatrix<double, 0, int> &);
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);

  double trace0(Eigen::SparseMatrix<double, 0, int> &);
  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return LU_K.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { return LU_K.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd v) { return LU_K.solve(v); }
  double logdet();
  Eigen::VectorXd rMVN(Eigen::VectorXd &, Eigen::VectorXd &)
  {
    std::cout << "lu_sparse_solver:rMVN not implimented\n";
    throw;
  };
};

#endif