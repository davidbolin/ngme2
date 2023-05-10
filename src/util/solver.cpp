#include "../include/solver.h"
using namespace Eigen;
// double myround(double x)
// {
//   if (x > 0)
//   {
//     return -1.0;
//   }
//   else
//   {
//     return 1.0;
//   }
// }

// void iterative_solver::initFromList(int nin, const Rcpp::List &init_list)
// {
//   if (nin <= 0)
//   {
//     std::cout << "iterative_solver::initFromList nin must be greater then >0\n";
//     throw("error");
//   }
//   n = nin;

//   N = Rcpp::as<int>(init_list["trace.iter"]);
//   R.setMaxIterations(Rcpp::as<int>(init_list["trace.solver.max.iter"]));
//   R.setTolerance(Rcpp::as<double>(init_list["trace.solver.tol"]));

//   // R.preconditioner().setDroptol(Rcpp::as<double>(init_list["solver.drop.tol"]));
//   // R.preconditioner().setFillfactor(Rcpp::as<double>(init_list["solver.fill.factor"]));

//   QU.setZero(n, N);
//   U.setRandom(n, N);
//   // U = U.unaryExpr(std::ptr_fun(myround));
//   U = U.unaryExpr([] (double x) -> double {myround(x);});

//   MQU.setZero(n, N);
// }

// void iterative_solver::init(int nin, int Nin, int max_iter, double tol)
// {
//   N = Nin;
//   n = nin;
//   R.setMaxIterations(max_iter);
//   R.setTolerance(tol);

//   QU.setZero(n, N);
//   U.setRandom(n, N);
//   // U = U.unaryExpr(std::ptr_fun(myround));
//   U = U.unaryExpr([] (double x) -> double {myround(x);});

//   MQU.setZero(n, N);
// }

// void iterative_solver::compute(SparseMatrix<double, 0, int> &M)
// {
//   R.factorize(M);
//   // U.setRandom(n,N);
//   // U = U.unaryExpr(std::ptr_fun(myround));
//   QU_computed = 0;
// }

// // Solve trace(M*Q^-1)
// double iterative_solver::trace(MatrixXd &M)
// {
//   if (QU_computed == 0)
//   {
//     for (int i = 0; i < N; i++)
//     {
//       QU.col(i) = R.solveWithGuess(U.col(i), QU.col(i));
//     }
//     QU_computed = 1;
//   }

//   MQU = M * QU;
//   double t = 0;

//   for (int i = 0; i < N; i++)
//   {
//     t += U.col(i).dot(MQU.col(i));
//   }

//   return t / N;
// }

// double iterative_solver::trace(SparseMatrix<double, 0, int> &M)
// {

//   if (QU_computed == 0)
//   {
//     for (int i = 0; i < N; i++)
//     {
//       QU.col(i) = R.solveWithGuess(U.col(i), QU.col(i));
//     }
//     QU_computed = 1;
//   }
//   MQU = M * QU;
//   double t = 0;
//   for (int i = 0; i < N; i++)
//   {
//     t += U.col(i).dot(MQU.col(i));
//   }
//   return t / N;
// }

// // sovle trace(M1*Q^-1*M2*Q^-1)
// double iterative_solver::trace2(SparseMatrix<double, 0, int> &M1, SparseMatrix<double, 0, int> &M2)
// {

//   if (QU_computed == 0)
//   {
//     for (int i = 0; i < N; i++)
//     {
//       QU.col(i) = R.solveWithGuess(U.col(i), QU.col(i));
//     }
//     QU_computed = 1;
//   }
//   MU = M2.transpose() * U;
//   MQU = M1 * QU;
//   for (int i = 0; i < N; i++)
//   {
//     MU.col(i) = R.solveWithGuess(MU.col(i), MU.col(i));
//   }
//   double t = 0;
//   for (int i = 0; i < N; i++)
//   {
//     t += MU.col(i).dot(MQU.col(i));
//   }
//   return t / N;
// }

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

// Sample Q^-1mu + chol(Q)Z = N(Q^-1mu,Q^-1)
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
  std::cout << "lu_solver::trace2 not implemented\n";
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
}

void lu_sparse_solver::initFromList(int nin, Rcpp::List const &)
{
  n = nin;
  KKtinv.resize(n, n);
  KKtinv_computed = 0;
}

void lu_sparse_solver::compute(const SparseMatrix<double, 0, int> &K_in)
{
  K = K_in;

  if (K.isCompressed() == 0)
    K.makeCompressed();

  if (K.rows() != n)
  {
    std::cout << "incorrect matrix size: n= " << n;
    std::cout << ", K = " << K.rows() << " * " << K.cols() << std::endl;
  }

  LU_K.factorize(K);
  L_KKt.factorize(K.transpose() * K); // KTK
  KKtinv_computed = 0;
}

// similar to compute
void lu_sparse_solver::computeKTK(const SparseMatrix<double, 0, int> &K_in)
{
  K = K_in;

  if (K.isCompressed() == 0)
    K.makeCompressed();

  if (K.rows() != n)
  {
    std::cout << "incorrect matrix size: n= " << n;
    std::cout << ", K = " << K.rows() << " * " << K.cols() << std::endl;
  }

  L_KKt.factorize(K.transpose() * K); // KTK
  KKtinv_computed = 0;
}

// Solve trace(K^-1 M)
double lu_sparse_solver::trace(const MatrixXd &M)
{

  std::cout << "lu_sparse_solver::trace for dense matrix M not implimented\n";
  throw;
  return 0;
}

double lu_sparse_solver::trace0(SparseMatrix<double, 0, int> &M)
{
  std::cout << "lu_sparse_solver0 not implemented";
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

double lu_sparse_solver::trace2(const SparseMatrix<double, 0, int> &M1, SparseMatrix<double, 0, int> &M2)
{
  std::cout << "lu_sparse_solver::trace2 not implimented\n";
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
}