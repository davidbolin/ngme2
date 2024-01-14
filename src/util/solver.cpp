#include "../include/solver.h"
using namespace Eigen;
double myround(double x){
	if(x > 0){
		return -1.0;
	} else {
		return 1.0;
	}
}

/* --------------------------------------------------------- */
// CHOLESKY SOLVER
void cholesky_solver::init(int nin, int Nin, int max_iter, double tol)
{
  n = nin;
  Qi.resize(n, n);
  Qi_computed = 0;
}

void cholesky_solver::initFromList(int nin, Rcpp::List const &)
{
  n = nin;
  Qi.resize(n, n);
  Qi_computed = 0;
}

void cholesky_solver::compute(const SparseMatrix<double, 0, int> &M)
{
  R.factorize(M);
  Qi_computed = 0;
  QU_computed = 0;
}

void cholesky_solver::set_ld()
{
  SparseMatrix<double, 0, int> R_Q = R.matrixL();
  ld = 2.0 * R_Q.diagonal().array().log().sum();
}

// Solve trace(M*Q^-1)
double cholesky_solver::trace(const MatrixXd &M)
{

  if (Qi_computed == 0)
  {
    SparseMatrix<double, 0, int> R_Q = R.matrixL();
    Qi = Qinv(R_Q);
    Qi_computed = 1;
  }

  MatrixXd Mreo;
  Mreo = R.permutationP().transpose() * M * R.permutationP();
  return Qi.cwiseProduct(Mreo).sum();
}

SparseMatrix<double, 0, int> cholesky_solver::return_Qinv()
{

  if (Qi_computed == 0)
  {
    SparseMatrix<double, 0, int> R_Q = R.matrixL();
    Qi = Qinv(R_Q);
    Qi_computed = 1;
  }

  SparseMatrix<double, 0, int> Qi_reo;
  Qi_reo = Qi.twistedBy(R.permutationPinv());
  return (Qi_reo);
}

double cholesky_solver::trace_num(const SparseMatrix<double, 0, int> &M)
{
  if (QU_computed==0){
    // U.setRandom(n,N);
    // The MatrixXd::Random() is complained by R CMD check
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < N; ++j) {
        U(i, j) = R::runif(-1.0, 1.0);
      }
    }

    U = U.unaryExpr(std::ref(myround));
    for(int i=0; i<N; i++){
      QU.col(i) = R.solve(U.col(i));
    }
    QU_computed = 1;
  }

  Eigen::MatrixXd MQU = M*QU;
  double t = 0;
  for(int i=0;i<N;i++){
    t += U.col(i).dot(MQU.col(i));
  }
  return t/N;
}

double cholesky_solver::trace(const SparseMatrix<double, 0, int> &M)
{
  if (Qi_computed == 0)
  {
    SparseMatrix<double, 0, int> R_Q = R.matrixL();
    Qi = Qinv(R_Q);
    Qi_computed = 1;
  }
  SparseMatrix<double, 0, int> Mreo;
  Mreo = M.twistedBy(R.permutationP());
  return Qi.cwiseProduct(Mreo).sum();
  /*
  double tr = 0.0;
  for(int i =0;i<n;i++){
    VectorXd r = M.col(i);
    VectorXd rv = R.solve(r);
    tr += rv[i];
  }
  return tr;
  */
}

double cholesky_solver::trace2(const SparseMatrix<double, 0, int> &M1, SparseMatrix<double, 0, int> &M2)
{
  /*
  std::cout << "trace2 is not working for cholesky!" << std::endl;
  if(Qi_computed == 0){
    SparseMatrix<double,0,int> R_Q = R.matrixL();
    Qi = Qinv(R_Q);
    Qi_computed = 1;
  }
  */
  SparseMatrix<double, 0, int> tmp1 = R.solve(M2);
  SparseMatrix<double, 0, int> tmp2 = R.solve(M1);
  tmp1 = tmp1 * tmp2;
  return tmp1.diagonal().sum();
}

// Sample Q^-1mu + chol(Q)Z = N(Q^-1mu,Q^-1), R is chol(Q)
Eigen::VectorXd cholesky_solver::rMVN(Eigen::VectorXd &mu, Eigen::VectorXd &z)
{

  // dest = R.permutationP() * mu;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> RP;
  RP.resize(mu.size());
  RP = R.permutationP();
  Eigen::VectorXd dest = R.permutationP() * mu;
  dest = R.matrixL().solve(dest);
  // here dest is the mean?
  dest = R.matrixU().solve(dest + z);
  dest = R.permutationPinv() * dest;
  return dest;
}

Eigen::VectorXd cholesky_solver::Qinv_diag()
{
  if (Qi_computed == 0)
  {
    SparseMatrix<double, 0, int> R_Q = R.matrixL();
    Qi = Qinv(R_Q);
    Qi_computed = 1;
  }
  VectorXd vars = Qi.diagonal();
  return R.permutationPinv() * vars;
}

/* --------------------------------------------------------- */
/* --------------------------------------------------------- */
// LU SOLVER
/* --------------------------------------------------------- */
/* --------------------------------------------------------- */
void lu_solver::init(int nin, int Nin, int max_iter, double tol)
{
  n = nin;
  Kinv.resize(n, n);
  Kinv_computed = 0;
}

void lu_solver::initFromList(int nin, Rcpp::List const &in_list)
{
  n = nin;
  Kinv.resize(n, n);
  Kinv_computed = 0;
}

void lu_solver::analyze(const Eigen::SparseMatrix<double, 0, int> &K)
{
}

void lu_solver::compute(const Eigen::SparseMatrix<double, 0, int> &M)
{
  LU_K.compute(MatrixXd(M));
}

void lu_solver::compute(const Eigen::MatrixXd &K)
{

  LU_K.compute(K);
  Kinv_computed = 0;
}

Eigen::VectorXd lu_solver::solve(Eigen::VectorXd &v, Eigen::VectorXd &d)
{
  return (LU_K.solve(v));
}

// Solve trace(M*K^-1)
double lu_solver::trace(const Eigen::MatrixXd &M)
{
  if (Kinv_computed == 0)
  {
    Kinv = LU_K.inverse();
    Kinv_computed = 1;
  }
  Eigen::MatrixXd tmp1 = M * Kinv;
  return (tmp1.diagonal().sum());
}

double lu_solver::trace(const Eigen::SparseMatrix<double, 0, int> &M)
{
  if (Kinv_computed == 0)
  {
    Kinv = LU_K.inverse();
    Kinv_computed = 1;
  }
  Eigen::MatrixXd tmp1 = M * Kinv;
  return (tmp1.diagonal().sum());
}

double lu_solver::trace2(const SparseMatrix<double, 0, int> &, SparseMatrix<double, 0, int> &g)
{
  // std::cout << "lu_solver::trace2 not implemented\n";
  throw;
}

inline double lu_solver::logdet()
{
  return (log(fabs(LU_K.determinant())));
}

/* --------------------------------------------------------- */
// LU sparse SOLVER

void lu_sparse_solver::init(int nin, int Nin, int max_iter, double tol)
{
  n = nin;
  KKtinv.resize(n, n);
  KKtinv_computed = 0;
  QU_computed = 0;
}

void lu_sparse_solver::initFromList(int nin, Rcpp::List const &)
{
  n = nin;
  KKtinv.resize(n, n);
  KKtinv_computed = 0;
  QU_computed = 0;
}

void lu_sparse_solver::compute(
  const SparseMatrix<double, 0, int> &K_in)
{
  K = K_in;

  if (K.isCompressed() == 0)
    K.makeCompressed();

  if (K.rows() != n)
  {
    // std::cout << "incorrect matrix size: n= " << n;
    // std::cout << ", K = " << K.rows() << " * " << K.cols() << std::endl;
  }

  LU_K.factorize(K);
  L_KKt.factorize(K.transpose() * K); // KTK
  KKtinv_computed = 0;
  QU_computed = 0;
}

// similar to compute
void lu_sparse_solver::compute_KTK(const SparseMatrix<double, 0, int> &K_in)
{
  K = K_in;

  if (K.isCompressed() == 0)
    K.makeCompressed();

  if (K.rows() != n)
  {
    // std::cout << "incorrect matrix size: n= " << n;
    // std::cout << ", K = " << K.rows() << " * " << K.cols() << std::endl;
  }

  L_KKt.factorize(K.transpose() * K);
  KKtinv_computed = 0;
}

void lu_sparse_solver::compute_LU(const SparseMatrix<double, 0, int> &K_in)
{
  K = K_in;

  if (K.isCompressed() == 0)
    K.makeCompressed();

  if (K.rows() != n)
  {
    // std::cout << "incorrect matrix size: n= " << n;
    // std::cout << ", K = " << K.rows() << " * " << K.cols() << std::endl;
  }

  LU_K.factorize(K);
  KKtinv_computed = 0;
  QU_computed = 0;
}

// Solve trace(K^-1 M)
double lu_sparse_solver::trace(const MatrixXd &M)
{

  // std::cout << "lu_sparse_solver::trace for dense matrix M not implimented\n";
  throw;
  return 0;
}

double lu_sparse_solver::trace0(SparseMatrix<double, 0, int> &M)
{
  // std::cout << "lu_sparse_solver0 not implemented";
  throw;
  return 0;
}

double lu_sparse_solver::trace(const SparseMatrix<double, 0, int> &M)
{
  if (KKtinv_computed == 0)
  {
    SparseMatrix<double, 0, int> L = L_KKt.matrixL();
    Qinv(L);
    KKtinv = Qinv(L);
    KKtinv_computed = 1;
  }
  SparseMatrix<double, 0, int> Mreo;
  Mreo = (K.transpose() * M).twistedBy(L_KKt.permutationP());
  return KKtinv.cwiseProduct(Mreo).sum();
}

double lu_sparse_solver::trace_num(const SparseMatrix<double, 0, int> &M)
{
  if (QU_computed==0) {
    // U.setRandom(n,N);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < N; ++j) {
        U(i, j) = R::runif(-1.0, 1.0);
      }
    }

    U = U.unaryExpr(std::ref(myround));
    for(int i=0; i<N; i++){
      QU.col(i) = LU_K.solve(U.col(i));
    }
    QU_computed = 1;
  }

  Eigen::MatrixXd MQU = M*QU;
  double t = 0;
  for(int i=0;i<N;i++){
    t += U.col(i).dot(MQU.col(i));
  }
  return t/N;
}

double lu_sparse_solver::trace2(const SparseMatrix<double, 0, int> &M1, SparseMatrix<double, 0, int> &M2)
{
  // std::cout << "lu_sparse_solver::trace2 not implimented\n";
  throw;
  return 0;
}

inline double lu_sparse_solver::logdet()
{
  SparseMatrix<double, 0, int> R_Q = L_KKt.matrixL();
  return R_Q.diagonal().array().log().sum();
}

void lu_sparse_solver::analyze(const Eigen::SparseMatrix<double, 0, int> &M)
{
  // if (M.isCompressed() == 0) M.makeCompressed();
  L_KKt.analyzePattern(M.transpose() * M);
  LU_K.analyzePattern(M);
  n=M.cols();
  U.resize(n, N);
  QU.resize(n, N);
  QU_computed = false;
}