#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <sstream>
#include "realtypes.h"

// three ways to initialize/ assign value to a matrix object:
// 1, matrix A; A.appendRow(std::vector<REAL>);	this is not a convenient way, 
// 2, matrix A(int i, int j); (=> Aij = 0) then we can use A(i,j)=REAL to assign value 
// each element of A
// 3, matrix A; A = ["1,2,2;2,3,4"] ([] is overloaded which will return a matrix), this is plan

// written on Feb 13, 2013
namespace dem{

class matrix{
      public:
	std::vector<std::vector<REAL> > value;
	int num_row;	// number of rows
	int num_col;	// number of columes

	matrix():num_row(0), num_col(0) {}	// constructor
	matrix(int, int);	// constructor, way 2
	std::vector<REAL> getCol(int i);	// get ith colume, 1<= i <=num_col
	std::vector<REAL> getRow(int i);	// get ith row
	void appendRow(std::vector<REAL>);
	void appendCol(std::vector<REAL>);
	void clear();	// clear all value of the matrix, assign num_row and num_col to be 0
	matrix getInvs();	// return the inverse 
	matrix getTrans();	// return the transpose
	REAL getNorm();	// return the norm of matrix, at present it can only calculate norm of a vector

	matrix& operator += (matrix A) 	{(*this) = (*this)+A; return *this;}
	matrix& operator += (REAL k) 	{(*this) = (*this)+k; return *this;}
	matrix& operator -= (matrix A) 	{(*this) = (*this)-A; return *this;}
	matrix& operator -= (REAL k) 	{(*this) = (*this)-k; return *this;}
	matrix& operator *= (matrix A) 	{(*this) = (*this)*A; return *this;}
	matrix& operator *= (REAL k)	{(*this) = (*this)*k; return *this;}
	matrix& operator = (matrix A);	// overload = for matrix
	REAL& operator () (int i, int j);	// access to the element at ith row and jth column: 1<= i <= num_rol
	std::string print();		// used for debug

	// non-member functions
	friend matrix operator + (matrix A, matrix B);	// overload + for two matrix
	friend matrix operator + (matrix, REAL);	// overload + for two matrix
	friend matrix operator + (REAL, matrix);	// overload + for scalar and matrix
	friend matrix operator - (matrix A, matrix B);	// overload - for matrix and scalar
	friend matrix operator - (REAL, matrix);	// overload - for scalar and matrix
	friend matrix operator - (matrix, REAL);	// overload - for matrix and scalar
	friend matrix operator * (matrix A, matrix B);	// overload * for two matrix
	friend matrix operator * (matrix A, REAL k);	// overload * for matrix and scalar
	friend matrix operator * (REAL k, matrix A);	// overload * for a scalar and a matrix
//	friend matrix operator / (matrix A, matrix B);
	friend matrix operator / (matrix A, REAL k);	// overload / for matrix and a scalar
//	friend matrix operator [] (stringstream);	// overload [] for assignment
	friend matrix expm(matrix);	// calculate the exponential of each element in the matrix. Not used here
	friend REAL trace(matrix);	// calculate the trace of the matrix, add on March 8, 2013 for fabric invariants

	

};
// all the operations defined in this matrix class can deal with any dimension matrix except getInvs(), which is now only
// suitable for 2 by 2 matrix and 1x1 matrix. In the future, I will expand it to any dimensions matrix using the method 
// taught in numerical method course named LU Decomposition


}//end of namespace dem

#endif
