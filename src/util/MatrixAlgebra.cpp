#include "../include/MatrixAlgebra.h"
#include <Rcpp.h>
#include <RcppEigen.h>

typedef Eigen::Triplet<double> T;
using namespace Eigen;
using namespace std;

MatrixXd communicationMatrix(const int n, const int m)
{
	MatrixXd K(n * m, n * m);

	K.setZero();
	for (int i = 0; i < n * m; i++)
		K(i, floor(i / m) + n * (i % m)) = 1;

	return K;
}

// Set Msub = M(ind,:) where ind is a 0-1 index vector
void get_submatrix(const MatrixXd &M, const VectorXi &ind, MatrixXd &Msub)
{
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "get_submatrix:: number of rows in M does not match length of index vector\n";
	}
	Msub.setZero(ind.sum(), M.cols());
	int k = 0;
	for (int i = 0; i < M.rows(); i++)
	{
		if (ind(i) == 1)
			Msub.row(k++) = M.row(i);
	}
}

// Return Msub = M(ind,:) where ind is a 0-1 index vector
MatrixXd get_submatrix(const MatrixXd &M, const VectorXi &ind)
{
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "get_submatrix:: number of rows in M does not match length of index vector\n";
	}
	MatrixXd Msub(ind.sum(), M.cols());
	int k = 0;
	for (int i = 0; i < M.rows(); i++)
	{
		if (ind(i) == 1)
			Msub.row(k++) = M.row(i);
	}
	return Msub;
}

MatrixXi get_submatrix(const MatrixXi &M, const VectorXi &ind)
{
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "get_submatrix:: number of rows in M does not match length of index vector\n";
	}
	MatrixXi Msub(ind.sum(), M.cols());
	int k = 0;
	for (int i = 0; i < M.rows(); i++)
	{
		if (ind(i) == 1)
		{
			Msub.row(k) = M.row(i);
			k++;
		}
	}
	return Msub;
}

void get_submatrix(SparseMatrix<double, 0, int> &M, const VectorXi &ind, SparseMatrix<double, 0, int> *Msub)
{
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "get_submatrix:: number of rows in M does not match length of index vector\n";
	}
	Msub->resize(ind.sum(), M.cols());
	// create vector with row indices for Msub
	VectorXi ind2(ind.size());
	int j = 0;
	for (int i = 0; i < ind.size(); i++)
	{
		j += ind(i);
		ind2(i) = j - 1.0;
	}
	for (int k = 0; k < M.outerSize(); ++k)
	{
		for (SparseMatrix<double, 0, int>::InnerIterator it(M, k); it; ++it)
		{
			if (ind(it.row()) == 1)
				Msub->insert(ind2(it.row()), it.col()) = it.value();
		}
	}
}

// set Ysub = Y(ind)
void get_subvector(const VectorXd &Y, const VectorXi &ind, VectorXd &Ysub)
{
	Ysub.resize(ind.sum());
	if (ind.size() != Y.size())
	{
		Rcpp::Rcout << "get_subvector:: number of rows in Y does not match length of index vector\n";
	}
	int k = 0;
	for (int i = 0; i < Y.size(); i++)
	{
		if (ind(i) == 1)
			Ysub(k++) = Y(i);
	}
}

// Set M(ind,:) = Msub where ind is a 0-1 index vector
void set_submatrix(MatrixXd &M, const MatrixXd &Msub, const VectorXi &ind)
{
	int k = 0;
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "set_submatrix:: number of rows in M does not match length of index vector\n";
	}
	for (int i = 0; i < ind.size(); i++)
	{
		if (ind(i) == 1)
		{
			M.row(i) = Msub.row(k);
			k++;
		}
	}
}

void set_submatrix(MatrixXi &M, const MatrixXi &Msub, const VectorXi &ind)
{
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "set_submatrix:: number of rows in M does not match length of index vector\n";
	}
	int k = 0;
	for (int i = 0; i < ind.size(); i++)
	{
		if (ind(i) == 1)
		{
			M.row(i) = Msub.row(k);
			k++;
		}
	}
}
// Set M(ind,r) = v
void set_subcol(MatrixXd &M, const int r, const VectorXi &ind, const VectorXd &v)
{
	int k = 0;
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "set_subcol:: number of rows in M does not match length of index vector\n";
	}
	for (int i = 0; i < ind.size(); i++)
	{
		if (ind(i) == 1)
		{
			M(i, r) = v(k);
			k++;
		}
	}
}

// Set M(ind,:) += Msub where ind is a 0-1 index vector
void add_submatrix(MatrixXd &M, const MatrixXd &Msub, const VectorXi &ind)
{
	int k = 0;
	if (ind.size() != M.rows())
	{
		Rcpp::Rcout << "add_submatrix:: number of rows in M does not match length of index vector\n";
	}
	for (int i = 0; i < ind.size(); i++)
	{
		if (ind(i) == 1)
		{
			M.row(i) += Msub.row(k);
			k++;
		}
	}
}

MatrixXi duplicatematrix(int n)
{
	MatrixXi D;
	D.setZero(n * n, n * (n - 1) / 2 + n);
	int k = n;
	int l = n;
	for (int i = 0; i < n; i++)
	{
		D(i * n + i, i) = 1;
		for (int j = i + 1; j < n; j++)
		{
			D(i * n + j, k++) = 1;
			D(j * n + i, l++) = 1;
		}
	}
	return D;
}

VectorXd vec(const MatrixXd &M)
{
	VectorXd V;
	V.setZero(M.rows() * M.cols());
	for (int i = 0; i < M.cols(); i++)
	{
		V.segment(i * M.rows(), M.rows()) = M.col(i);
	}
	return V;
}

MatrixXd veci(VectorXd &v, int n, int m)
{
	if (v.size() != n * m)
	{
		cout << "Wrong dimensions in reshape: " << v.size() << " (" << n << "," << m << ")" << endl;
	}
	MatrixXd M;
	M.setZero(n, m);
	int count = 0;
	for (int i = 0; i < m; i++)
	{
		M.col(i) = v.segment(count, n);
		count += n;
	}
	return M;
}

VectorXd vech(MatrixXd &M)
{
	int n = M.rows();
	VectorXd V;
	V.setZero(n * (n - 1) / 2 + n);
	V.segment(0, n) = M.diagonal();
	int k = n;
	for (int i = 0; i < n - 1; i++)
	{
		V.segment(k, n - i - 1) = M.block(i, i + 1, 1, n - i - 1).transpose();
		k += n - i - 1;
	}
	return V;
}

void read_SparseMatrix(SparseMatrix<double, 0, int> &M, string path, string name)
{

	FILE *pFile;
	pFile = fopen((path + name + "i.bin").c_str(), "rb");
	if (pFile == NULL)
	{
		fputs("File error", stderr);
		exit(1);
	}

	int n, m, ni;
	if (1 != fread(&n, sizeof(int), 1, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	if (1 != fread(&m, sizeof(int), 1, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}

	if (1 != fread(&ni, sizeof(int), 1, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	M.resize(n, m);
	MatrixXi ij;
	ij.setZero(ni, 2);
	if (ni * 2 != fread(ij.data(), sizeof(int), ni * 2, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	fclose(pFile);
	VectorXd v;
	v.setZero(ni);
	pFile = fopen((path + name + "v.bin").c_str(), "rb");
	if (pFile == NULL)
	{
		fputs("File error", stderr);
		exit(1);
	}
	if (ni != fread(v.data(), sizeof(double), ni, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	fclose(pFile);
	vector<T> coef;
	for (int i = 0; i < ni; i++)
	{
		coef.push_back(T(ij(i, 0), ij(i, 1), v(i)));
	}
	M.setFromTriplets(coef.begin(), coef.end());
}

void read_MatrixXd(MatrixXd &M, string path, string name)
{
	FILE *pFile;
	pFile = fopen((path + name + ".bin").c_str(), "rb");
	if (pFile == NULL)
	{
		fputs("File error", stderr);
		exit(1);
	}

	int n, d;
	if (1 != fread(&n, sizeof(int), 1, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	if (1 != fread(&d, sizeof(int), 1, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	M.setZero(n, d);
	if (n * d != fread(M.data(), sizeof(double), n * d, pFile))
	{
		fputs("Read error\n", stderr);
		exit(1);
	}
	fclose(pFile);
}

SparseMatrix<double, 0, int> Qinv(SparseMatrix<double, 0, int> &Q_R)
{

	int *Rir;
	int *Rjc;
	double *Rpr;
	Rir = Q_R.innerIndexPtr();
	Rjc = Q_R.outerIndexPtr();
	Rpr = Q_R.valuePtr();
	int n = Q_R.rows();
	typedef pair<size_t, double> Qpairtype;
	typedef vector<Qpairtype> Qvectype;
	typedef vector<Qvectype> Qtype;

	int i, j;

	/*
	 * Copy cholesky factor to more convenient format and extract diagonal elements
	 */

	Qtype R(n);
	vector<double> D(n);
	// Extract the elements and store the sparse R-matrix
	// in a more convinient format.
	if (Rjc[n] - Rjc[n - 1] == 1)
	{
		// only one element in the last column, assume lower triangular matrix
		for (int c = 0; c < n; ++c)
		{
			D[c] = Rpr[Rjc[c]];
			R[c].resize(Rjc[c + 1] - Rjc[c]);
			for (j = Rjc[c], i = 0; j < Rjc[c + 1]; ++j, ++i)
				R[c][i] = Qpairtype(Rir[j], Rpr[j]);
		}
	}
	else
	{
		// assume upper triangular matrix - first find number of element in each row
		vector<size_t> nRow(n), iRow(n);
		for (int c = 0; c < n; ++c)
		{
			for (j = Rjc[c]; j < Rjc[c + 1]; ++j)
				++nRow[Rir[j]];
			D[c] = Rpr[Rjc[c + 1] - 1];
		}
		for (int c = 0; c < n; ++c)
			R[c].resize(nRow[c]);
		for (int c = 0; c < n; ++c)
		{
			for (j = Rjc[c]; j < Rjc[c + 1]; ++j)
				R[Rir[j]][iRow[Rir[j]]++] = Qpairtype(c, Rpr[j]);
		}
	}

	/*
	 * Calculate inverse
	 */
	Qvectype::iterator pos;
	size_t Nmax = 0;
	// divide all elemnts in R by the diagonal-elements
	for (i = 0; i < n; ++i)
	{
		// find the maximal number of non-zero elements in any row of R
		if (Nmax < R[i].size())
			Nmax = R[i].size();
		// compute R[i,j]/D[i]
		for (pos = R[i].begin(); pos != R[i].end(); ++pos)
			(pos->second) /= D[i];
		// and compute 1/d^2
		D[i] = 1 / (D[i] * D[i]);
	}

	// count number of elements that is going to end up in iQ
	vector<size_t> nnz(n, 1);
	for (i = 0; i < n; ++i)
	{
		// first find the indices of the non-zero elements
		for (pos = R[i].begin(), ++pos; pos != R[i].end(); ++pos)
		{
			nnz[i]++;
			nnz[pos->first]++;
		}
	}

	// vectors containing the location and values within one column
	vector<size_t> ii(Nmax);
	vector<double> s(Nmax);
	vector<Qvectype::iterator> iQpos(Nmax);
	vector<Qvectype::iterator> iQstart(n);

	// create a structure holding the inverse matrix
	Qtype iQ(n);
	for (i = 0; i < n; ++i)
	{
		iQ[i].resize(nnz[i]);
		iQstart[i] = iQ[i].end();
	}
	// loop over the columns of the matrix
	i = n;
	while (i > 0)
	{
		--i;
		// first find the indices of the non-zero elements
		for (pos = R[i].begin(), ++pos, j = 0; pos != R[i].end(); ++pos, j++)
		{
			ii[j] = pos->first;		   // index of elements
			s[j] = 0;				   // set values to zero
			iQpos[j] = iQstart[ii[j]]; // start of each iQ row
		}
// multiply the row of R with the rows of iQ
#pragma omp parallel for private(pos)
		for (int j2 = 0; j2 < (R[i].size() - 1); ++j2)
		{
			Qvectype::iterator iQpos_tmp = iQpos[j2];
			Qvectype::iterator iQend = iQ[ii[j2]].end();
			for (pos = R[i].begin(), ++pos; pos != R[i].end(); ++pos)
			{
				for (; iQpos_tmp != iQend && iQpos_tmp->first < pos->first; ++iQpos_tmp)
				{
				}
				if (iQpos_tmp != iQend && iQpos_tmp->first == pos->first)
					s[j2] += (iQpos_tmp->second) * (pos->second);
			}
		}
		// the diagonal elements
		double diag = D[i];
		for (pos = R[i].begin(), ++pos, j = 0; pos != R[i].end(); ++pos, ++j)
			diag += s[j] * (pos->second);

		// add the elements to iQ
		j = R[i].size() - 1;
		while (j > 0)
		{
			--j;
			*(--iQstart[i]) = Qpairtype(ii[j], -s[j]);
			*(--iQstart[ii[j]]) = Qpairtype(i, -s[j]);
		}
		*(--iQstart[i]) = Qpairtype(i, diag);
	}
	// construct sparse matrix.
	SparseMatrix<double, 0, int> Qinverse;
	Qinverse.resize(n, n);
	vector<T> coef;

	for (i = 0, j = 0; i < n; ++i)
	{
		for (pos = iQ[i].begin(); pos != iQ[i].end(); ++pos, ++j)
		{
			coef.push_back(T(i, pos->first, pos->second));
		}
	}
	Qinverse.setFromTriplets(coef.begin(), coef.end());
	return Qinverse;
}

VectorXi ind2sub(int k, int nrows, int ncols)
{
	VectorXi v;
	v.setZero(2);
	v(1) = floor(k / nrows);
	v(0) = k - v(1) * nrows;
	return v;
}

SparseMatrix<double, 0, int> Qinv2(SparseMatrix<double, 0, int> &L)
{
	int n = L.rows();
	SparseMatrix<double, 0, int> Sigma;
	Sigma.resize(n, n);
	Sigma.reserve(L.nonZeros() * 50);
	int j;
	double Lii, val;
	for (int i = n - 1; i >= 0; i--)
	{
		for (SparseMatrix<double>::ReverseInnerIterator ij(L, i); ij; --ij)
		{

			j = ij.row();
			val = 0;
			SparseMatrix<double>::InnerIterator iL(L, i);
			SparseMatrix<double>::InnerIterator iS(Sigma, j);
			for (; iL; ++iL)
			{
				if (iL.row() == iL.col())
				{
					Lii = iL.value();
				}
				else
				{
					for (; iS.row() < iL.row(); ++iS)
					{
					}
					if (iS.row() == iL.row())
						val += iL.value() * iS.value();
				}
			}
			if (i == j)
			{
				Sigma.insert(i, j) = 1 / (Lii * Lii) - val / Lii;
			}
			else
			{
				Sigma.insert(i, j) = -val / Lii;
				Sigma.insert(j, i) = -val / Lii;
			}
		}
	}
	return Sigma;
}

SparseMatrix<double, 0, int> kronecker(SparseMatrix<double, 0, int> &A, SparseMatrix<double, 0, int> &B)
{
	int Br = B.rows();
	int Bc = B.cols();
	int Ar = A.rows();
	int Ac = A.cols();

	SparseMatrix<double, 0, int> AB;
	AB.resize(Ar * Br, Ac * Bc);

	VectorXi nnzA = VectorXi::Zero(Ac);
	for (int kA = 0; kA < A.outerSize(); ++kA)
	{
		for (SparseMatrix<double, 0, int>::InnerIterator itA(A, kA); itA; ++itA)
		{
			nnzA(itA.col())++;
		}
	}
	VectorXi nnzB = VectorXi::Zero(Bc);
	for (int kB = 0; kB < B.outerSize(); ++kB)
	{
		for (SparseMatrix<double, 0, int>::InnerIterator itB(B, kB); itB; ++itB)
		{
			nnzB(itB.col())++;
		}
	}
	Matrix<int, Dynamic, Dynamic, ColMajor> nnzAB = nnzB * nnzA.transpose();
	AB.reserve(VectorXi::Map(nnzAB.data(), nnzAB.size()));

	for (int kA = 0; kA < A.outerSize(); ++kA)
	{
		for (int kB = 0; kB < B.outerSize(); ++kB)
		{
			for (SparseMatrix<double, 0, int>::InnerIterator itA(A, kA); itA; ++itA)
			{
				for (SparseMatrix<double, 0, int>::InnerIterator itB(B, kB); itB; ++itB)
				{
					int i = itA.row() * Br + itB.row();
					int j = itA.col() * Bc + itB.col();
					AB.insert(i, j) = itA.value() * itB.value();
				}
			}
		}
	}
	return AB;
}

/*
  B of size n1 x n2
 Set A(i:i+n1,j:j+n2) = B
*/
void setSparseBlock(SparseMatrix<double, 0, int> *A, int i, int j, SparseMatrix<double, 0, int> &B)
{
	for (int k = 0; k < B.outerSize(); ++k)
	{
		for (SparseMatrix<double, 0, int>::InnerIterator it(B, k); it; ++it)
		{
			A->insert(it.row() + i, it.col() + j) = it.value();
		}
	}
}
/*
  B of size n1 x n2
 Set A(i:i+n1,j:j+n2) = B (update)
*/
void setSparseBlock_update(SparseMatrix<double, 0, int> *A, int i, int j, SparseMatrix<double, 0, int> &B)
{
	for (int k = 0; k < B.outerSize(); ++k)
	{
		for (SparseMatrix<double, 0, int>::InnerIterator it(B, k); it; ++it)
		{
			A->coeffRef(it.row() + i, it.col() + j) = it.value();
		}
	}
}

SparseMatrix<double, 0, int> full2sparse(MatrixXd &M)
{
	SparseMatrix<double, 0, int> A;
	A.resize(M.rows(), M.cols());
	for (int i = 0; i < M.rows(); i++)
	{
		for (int j = 0; j < M.cols(); j++)
		{
			A.insert(i, j) = M(i, j);
		}
	}
	return A;
}

Eigen::MatrixXd NormalOuterExpectation(const Eigen::MatrixXd &Sigma,
									   const Eigen::VectorXd &mu,
									   const Eigen::VectorXd &m)
{

	Eigen::MatrixXd Res;
	Res.setZero(Sigma.rows(), Sigma.cols());
	Eigen::VectorXd mu_m = mu - m;
	Res += Sigma;
	Res += mu_m * mu_m.transpose();
	return (Res);
}

Eigen::MatrixXd NormalFourthExpectation(const Eigen::MatrixXd &Sigma,
										const Eigen::VectorXd &mu,
										const Eigen::VectorXd &m)
{
	Eigen::MatrixXd Res;
	Res.setZero(Sigma.rows(), Sigma.cols());
	Eigen::MatrixXd mu_m_outer = (mu + m) * (mu + m).transpose();
	Eigen::MatrixXd Temp = (Sigma + mu_m_outer);
	Res += Temp * Temp;
	Res.array() *= 2 + Sigma.diagonal().sum();
	Res += mu_m_outer * (Sigma - mu_m_outer);
	return (Res);
}

Eigen::MatrixXd NormalouterKron(const Eigen::MatrixXd &Sigma,
								const Eigen::VectorXd &mu)
{
	Eigen::MatrixXd Res;
	int n = mu.size();
	Eigen::MatrixXd I;
	I.setIdentity(n * n, n * n);
	Eigen::MatrixXd Kn;
	Kn = communicationMatrix(n, n);

	Res = kroneckerProduct(Sigma, Sigma).eval();
	Eigen::MatrixXd mu_mu = mu * mu.transpose();
	Res += kroneckerProduct(Sigma, mu_mu).eval();
	Res += kroneckerProduct(mu_mu, Sigma).eval();
	Eigen::MatrixXd VX = (I + Kn) * Res;
	Eigen::MatrixXd Exxt = Sigma + mu_mu;
	Eigen::VectorXd EX = vec(Exxt);
	return (VX + EX * EX.transpose());
}

Eigen::MatrixXd Normalxxty(const Eigen::MatrixXd &Sigma1,
						   const Eigen::MatrixXd &Sigma12,
						   const Eigen::VectorXd &mu1,
						   const Eigen::VectorXd &mu2)
{

	Eigen::MatrixXd mu1Kmu1 = kroneckerProduct(mu1, mu1);
	Eigen::MatrixXd mu1Kmu1mu2t = kroneckerProduct(mu2.transpose(), mu1Kmu1);
	Eigen::MatrixXd Res = mu1Kmu1mu2t;
	Res += kroneckerProduct(mu1, Sigma12) + kroneckerProduct(Sigma12, mu1);
	Eigen::VectorXd vSigma1 = vec(Sigma1);
	Res += vSigma1 * mu2.transpose();
	return (Res);
}

Eigen::MatrixXd HessianSigma(const Eigen::MatrixXd &Z,
							 const Eigen::MatrixXd &iSigma,
							 const Eigen::MatrixXd &Sigma,
							 const int n)
{
	int d = Sigma.rows();
	Eigen::MatrixXd Kdd = communicationMatrix(d, d);
	Kdd.array() *= 0.5;
	Eigen::MatrixXd QZQ = iSigma * Z * iSigma;
	Eigen::MatrixXd Res = -kroneckerProduct(QZQ, iSigma);
	QZQ -= n * iSigma;
	Res -= kroneckerProduct(iSigma, QZQ);

	return (Kdd * Res);
}
