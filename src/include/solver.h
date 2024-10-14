#ifndef __Solver__Solver__
#define __Solver__Solver__

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#undef COMPLEX

#include <iostream>
#include <string.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include "MatrixAlgebra.h"


class solver
{
private:
public:
  int n;
  solver(const solver &){};
  solver(){};
  virtual ~solver(){};
  virtual void init(int, int, int, double) = 0;
  // virtual void initFromList(int, Rcpp::List const &) = 0;
  virtual inline void analyze(const Eigen::SparseMatrix<double, 0, int> &) = 0;
  virtual void compute(const Eigen::SparseMatrix<double, 0, int> &) = 0;
  virtual void compute(const Eigen::MatrixXd &) { std::cout << "compute not implimented for MatrixXd\n"; };
  virtual inline Eigen::VectorXd solve(Eigen::VectorXd &v) = 0;
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
  virtual Eigen::VectorXd rMVN(Eigen::VectorXd &, Eigen::VectorXd &) {
    std::cout << "rMVN not implimented\n";
    throw;
  }
  
  // sample from N(Q^-1 mu, Q^-1), Q = G^T G + H^T H
  virtual Eigen::VectorXd rMVN(
    SparseMatrix<double, 0, int> & G,
    SparseMatrix<double, 0, int> & H,
    Eigen::VectorXd & mu, 
    Eigen::VectorXd & z1, 
    Eigen::VectorXd & z2
  ) {
    Eigen::VectorXd x = G.transpose() * z1 + H.transpose() * z2 + mu;
    return solve(x);
  }

  virtual SparseMatrix<double, 0, int> return_Qinv()
  {
    SparseMatrix<double, 0, int> a;
    std::cout << "return_Qinv not implimented for MatrixXd\n";
    return a;
  }
};


class cholesky_solver : public virtual solver
{
private:
  bool Qi_computed {false}, QU_computed {false};
  double ld;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int> > R;
  Eigen::SparseMatrix<double, 0, int> Qi;
  Eigen::MatrixXd U, QU;
  int N {10};
  void set_ld();

public:
  cholesky_solver(const cholesky_solver &){};
  cholesky_solver(){};
  ~cholesky_solver(){};
  void init(int, int, int, double);
  // void initFromList(int, Rcpp::List const &);

  inline void set_N(int n) { N = n; }
  inline void analyze(const Eigen::SparseMatrix<double, 0, int> &M) {
    R.analyzePattern(M);
    QU_computed = false;
  }
  void compute(const Eigen::SparseMatrix<double, 0, int> &);

  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { 
    return R.solve(v); 
  }

  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { 
    return R.solve(v); 
  }

  inline Eigen::SparseMatrix<double, 0, int> solveMatrix(const Eigen::SparseMatrix<double, 0, int> &v) { return R.solve(v); }
  double trace(const Eigen::MatrixXd &);
  double trace(const Eigen::SparseMatrix<double, 0, int> &);
  double trace_num(const Eigen::SparseMatrix<double, 0, int> &);
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);
  inline double logdet()
  {
    set_ld();
    return ld;
  }
  VectorXd Qinv_diag();
  Eigen::VectorXd rMVN(Eigen::VectorXd &, Eigen::VectorXd &);
  Eigen::VectorXd rMVN(
    SparseMatrix<double, 0, int> & G,
    SparseMatrix<double, 0, int> & H,
    Eigen::VectorXd & mu, 
    Eigen::VectorXd & z1, 
    Eigen::VectorXd & z2
  );

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
};

class lu_sparse_solver : public virtual solver
{
private:
  int n;
  Eigen::SparseMatrix<double, 0, int> KKtinv;
  Eigen::SparseMatrix<double, 0, int> K;
  Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int> > LU_K;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int> > L_KKt;
  int KKtinv_computed, QU_computed;
  Eigen::MatrixXd U, QU;
  int N {10};

public:
  inline void set_N(int n) { N = n; }
  void init(int, int, int, double);
  // void initFromList(int, Rcpp::List const &);
  void analyze(const Eigen::SparseMatrix<double, 0, int> &);
  void compute(const Eigen::SparseMatrix<double, 0, int> &);
  void compute_KTK(const Eigen::SparseMatrix<double, 0, int> &);
  void compute_LU(const Eigen::SparseMatrix<double, 0, int> &);
  double trace(const Eigen::MatrixXd &);
  double trace(const Eigen::SparseMatrix<double, 0, int> &);
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);
  double trace_num(const Eigen::SparseMatrix<double, 0, int> &);

  double trace0(Eigen::SparseMatrix<double, 0, int> &);
  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return LU_K.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { return LU_K.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd v) { return LU_K.solve(v); }
  double logdet();
};


class iterative_solver : public virtual solver
{
private:
  int N;
  Eigen::MatrixXd U, QU, MQU, MU, prevQU;
  bool QU_computed;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double, 0, int>, Lower, Eigen::IncompleteCholesky<double> > R;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double, 0, int>, Lower|Upper, DiagonalPreconditioner<double> > R;
  int curr_iter {0};
public:
  ~iterative_solver(){};
  void init(int, int, int, double);
  void initFromList(int, Rcpp::List const &);
  
  // Update this line to match the base class signature
  inline void analyze(const Eigen::SparseMatrix<double, 0, int> &M) override { 
    R.compute(M);
    QU_computed = false;
  }
  
  void compute(const Eigen::SparseMatrix<double, 0, int> &);
  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { return R.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return R.solveWithGuess(v, x); }
  double trace(const Eigen::MatrixXd &);
  double trace(const Eigen::SparseMatrix<double, 0, int> &);
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &);
  double trace_num(const SparseMatrix<double, 0, int> &);
};

// if system is macOS, use Eigen::AccelerateLLT
#ifdef __APPLE__
#include <Eigen/AccelerateSupport>
class accel_llt_solver : public virtual solver
{
private:
  Eigen::AccelerateLLT<Eigen::SparseMatrix<double, 0, int> > R;
  bool Qi_computed {false}, QU_computed {false}; 
  Eigen::SparseMatrix<double, 0, int> Qi;
  Eigen::MatrixXd U, QU;
  int N {10};
public:
  void init(int, int, int, double);
  void analyze(const Eigen::SparseMatrix<double, 0, int> & M) {
    R.analyzePattern(M);
    QU_computed = false;
  }
  void compute(const Eigen::SparseMatrix<double, 0, int> & M) {
    R.factorize(M);
    Qi_computed = false;
    QU_computed = false;
  }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return R.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { return R.solve(v); }
  inline Eigen::SparseMatrix<double, 0, int> solveMatrix(const Eigen::SparseMatrix<double, 0, int> &v) { return R.solve(v); }
  
  double trace(const Eigen::MatrixXd &) {return 0.0;}
  double trace(const Eigen::SparseMatrix<double, 0, int> &) {return 0.0;}
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &) {return 0.0;}
  double trace_num(const SparseMatrix<double, 0, int> &);
};
#endif


#ifdef USEMKL
#include <Eigen/PardisoSupport>
class pardiso_llt_solver : public virtual solver {
private:
  Eigen::PardisoLLT<Eigen::SparseMatrix<double, 0, int> > R;
  bool Qi_computed {false}, QU_computed {false}; 
  Eigen::SparseMatrix<double, 0, int> Qi;
  Eigen::MatrixXd U, QU;
  int N {10};
public:
  void init(int, int, int, double);
  void analyze(const Eigen::SparseMatrix<double, 0, int> & M) {
    R.analyzePattern(M);
    QU_computed = false;
  }
  void compute(const Eigen::SparseMatrix<double, 0, int> & M) {
    R.factorize(M);
    Qi_computed = false;
    QU_computed = false;
  }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v, Eigen::VectorXd &x) { return R.solve(v); }
  inline Eigen::VectorXd solve(Eigen::VectorXd &v) { return R.solve(v); }
  inline Eigen::SparseMatrix<double, 0, int> solveMatrix(const Eigen::SparseMatrix<double, 0, int> &v) { return R.solve(v); }
  double trace(const Eigen::MatrixXd &) {return 0.0;}
  double trace(const Eigen::SparseMatrix<double, 0, int> &) {return 0.0;}
  double trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &) {return 0.0;}
  double trace_num(const SparseMatrix<double, 0, int> &);
};
#endif

#endif