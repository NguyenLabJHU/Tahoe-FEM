#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <sstream>
#include <math.h>
#include "realtypes.h"

// three ways to initialize/ assign value to a matrix object:
// 1, matrix A; A.appendRow(std::vector<REAL>);	this is not a convenient way, 
// 2, matrix A(int i, int j); (=> Aij = 0) then we can use A(i,j)=REAL to assign value 
// each element of A
// 3, matrix A; A = ["1,2,2;2,3,4"] ([] is overloaded which will return a matrix), this is plan

namespace memFluid{

class matrix{
      public:
	std::vector<std::vector<REAL> > value;	// matrix is vector of row vector
	int num_row;	// number of rows
	int num_col;	// number of columes

	matrix():num_row(0), num_col(0) {}	// constructor
	matrix(int, int);	// constructor, way 2
	void calcDimensions();
	std::vector<REAL> getCol(int i) const;	// get ith colume, 1<= i <=num_col
	std::vector<REAL> getRow(int i) const;	// get ith row
	void appendRow(const std::vector<REAL> &);
	void appendCol(const std::vector<REAL> &);
	void clear();	// clear all value of the matrix, assign num_row and num_col to be 0
	matrix getInvs();	// return the inverse 
	matrix getTrans() const;	// return the transpose
	REAL getNorm();	// return the norm of matrix, at present it can only calculate norm of a vector
        //void add(REAL);	// row vector
	void LU(matrix &L, matrix &U) const;	// LU decomposition, November 24
	REAL getVal(int , int) const;

	matrix& operator += (const matrix & A) 	{(*this) = (*this)+A; return *this;}
	matrix& operator += (REAL k) 	{(*this) = (*this)+k; return *this;}
	matrix& operator -= (const matrix & A) 	{(*this) = (*this)-A; return *this;}
	matrix& operator -= (REAL k) 	{(*this) = (*this)-k; return *this;}
	matrix& operator *= (const matrix & A) 	{(*this) = (*this)*A; return *this;}
	matrix& operator *= (REAL k)	{(*this) = (*this)*k; return *this;}
	matrix& operator = (const matrix &A);	// overload = for matrix
	REAL& operator () (int i, int j);	// access to the element at ith row and jth column: 1<= i <= num_rol
	std::string print() const;		// used for debug


// non-member functions
friend matrix operator + (const matrix &A, const matrix &B);	// overload + for two matrix
friend matrix operator + (const matrix &, REAL);	// overload + for two matrix
friend matrix operator + (REAL, const matrix &);	// overload + for scalar and matrix
friend matrix operator - (const matrix &A, const matrix &B);	// overload - for matrix and scalar
friend matrix operator - (REAL, const matrix &);	// overload - for scalar and matrix
friend matrix operator - (const matrix &, REAL);	// overload - for matrix and scalar
friend matrix operator * (const matrix &A, const matrix &B);	// overload * for two matrix
friend matrix operator * (const matrix &A, REAL k);	// overload * for matrix and scalar
friend matrix operator * (REAL k, const matrix &A);	// overload * for a scalar and a matrix
//friend matrix operator / (matrix A, matrix B);
friend matrix operator / (const matrix &A, REAL k);	// overload / for matrix and a scalar
friend matrix operator % (matrix &A, matrix &B);	// November 15, leftdivision, "\"

//friend matrix operator [] (stringstream);	// overload [] for assignment
friend matrix expm(const matrix &);	// calculate the exponential of each element in the matrix. Not used here

}; // end matrix

// October 31, 2013


bool isnan_mat(matrix &);	
int size(const matrix &, int);	// the same size function in matlab
matrix ones(int, int);		// the same ones function in matlab
matrix zeros(int, int);		// the same zeros function in matlab
REAL max(matrix &);	// the same max function in matlab
REAL min(matrix &);	// the same min function in matlab
matrix abs(matrix &);	// the same abs function in matlab
int length(const matrix &);	// the same length function in matlab
REAL norm(matrix &);	// the same norm function in matlab
REAL det(matrix &);	// the same det function in matlab
matrix linspace(REAL, REAL, int num);

void matrixEqnSolver(matrix &dsol, matrix &K, matrix &F);

}//end of namespace dem

#endif
