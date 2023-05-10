#ifndef __MatrixAlgebra__MatrixAlgebra__
#define __MatrixAlgebra__MatrixAlgebra__


#include <iostream>
#include <string.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/SparseCore>
#ifdef _OPENMP
	#include<omp.h>
#endif
using namespace Eigen;

/*
	Creating the communcation matrix,
	The matrix Knm s.t for an mxn A matrix
	Knm vec(A) =vec(A')

	output:
		Knm		(n*m)x(n*m( communcation matrix )

	reference : Magnus, Matrix Algebra p.331


 */
MatrixXd communicationMatrix(const int , const int );


/*
	Extract submatrix M(ind,:)
*/
MatrixXd get_submatrix(const MatrixXd&, const VectorXi&);
void get_submatrix(const MatrixXd&, const VectorXi&,MatrixXd&);
MatrixXi get_submatrix(const MatrixXi&, const VectorXi&);
void get_submatrix(SparseMatrix<double,0,int>&,const VectorXi&, SparseMatrix<double,0,int>*);

/*
	Set M(ind,:) = Msub
*/
void set_submatrix(MatrixXd&, const MatrixXd&, const VectorXi&);
void set_submatrix(MatrixXi&, const MatrixXi&, const VectorXi&);

/*
	M(ind,:) += Msub
*/
void add_submatrix(MatrixXd& M, const MatrixXd& Msub, const VectorXi& ind);

//set Ysub = Y(ind)
void get_subvector(const VectorXd&, const VectorXi&, VectorXd&);

// Set M(ind,r) = v
void set_subcol(MatrixXd&, const int, const VectorXi&, const VectorXd&);


// transform vech(M) to vec(M)
MatrixXd duplicatematrix(int);

//
VectorXd vec(const MatrixXd&);
VectorXd vech(MatrixXd&);
MatrixXd veci(VectorXd&,int,int);

SparseMatrix<double,0,int> Qinv(SparseMatrix<double,0,int>&);

VectorXi ind2sub(int,int,int);

void read_SparseMatrix(SparseMatrix<double,0,int>&,std::string,std::string);
void read_MatrixXd(MatrixXd&,std::string,std::string);

SparseMatrix<double,0,int> Qinv2(SparseMatrix<double,0,int>& L);

SparseMatrix<double,0,int> kronecker(SparseMatrix<double,0,int>&,SparseMatrix<double,0,int>&);
void setSparseBlock(SparseMatrix<double,0,int>*,int, int, const SparseMatrix<double,0,int>&);
void setSparseBlock_update(SparseMatrix<double,0,int>*,int, int, const SparseMatrix<double,0,int>&);

//convert a full matrix to sparse format
SparseMatrix<double,0,int> full2sparse(MatrixXd&);

// Computing the expectation of
// (x-m)(x-m)' when x ~  N(mu, Sigma)
Eigen::MatrixXd NormalOuterExpectation(const Eigen::MatrixXd &,
									   const Eigen::VectorXd &,
									   const Eigen::VectorXd &);
// Computing the expectation of
// xx'xx' when x ~  N(mu, Sigma)
Eigen::MatrixXd NormalFourthExpectation(const Eigen::MatrixXd &,
									    const Eigen::VectorXd &,
									    const Eigen::VectorXd &);
// Computing the expectation of
// xx' \otimes xx' when x ~  N(mu, Sigma)
Eigen::MatrixXd NormalouterKron(const Eigen::MatrixXd & ,
							    const Eigen::VectorXd  & );

// Computing the expectation of
// xx' \otimes xx' when [x,y] ~  N([mu_1, \mu_2], [Sigma_1, \Sigma_{12}, \Sigma_{21}, \Sigma_2])
Eigen::MatrixXd Normalxxty(const Eigen::MatrixXd & ,
						   const Eigen::MatrixXd & ,
						   const Eigen::VectorXd & ,
						   const Eigen::VectorXd & );
//
// computing the hessian of multivariate normal (not elementwise)
// Z      - XX'
// iSigma -  Sigma^{-1}
//  Sigma -  Sigma
Eigen::MatrixXd HessianSigma(const Eigen::MatrixXd & ,
						   const Eigen::MatrixXd & ,
						   const Eigen::MatrixXd & ,
						   const int );

SparseMatrix<double,0,int> kroneckerEigen(SparseMatrix<double,0,int>&,SparseMatrix<double,0,int>&);

#endif /* defined(__MatrixAlgebra__MatrixAlgebra__) */