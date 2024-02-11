#include "matrix.h"
#include <cstdlib>
#include <cmath>
// written on Feb 13, 2013
namespace dem{

matrix::matrix(int i, int j){
	if(i<1 || j<1){
		std::cout << "Error: the dimensions of matrix should be positive!" << std::endl;
		exit(-1);	
	}
	else{
		num_row = 0;
		num_col = 0;	// to avoid errors in appendRow()
		std::vector<REAL> temp_vec;
		for(int ir=0; ir!=i; ir++){
			temp_vec.clear();	// begin a new row			
			for(int ic=0; ic!=j; ic++){
				temp_vec.push_back(0);
			}
			appendRow(temp_vec);
		}
		num_row = i;
		num_col = j;
	}
}

std::vector<REAL> matrix::getCol(int i){
	if(i<=0 || i>num_col){
		std::cout << "Error: index exceeds!"	<< std::endl;
		exit(-1);
	}
	std::vector<REAL> result;
	result.clear();
	for(std::vector<std::vector<REAL> >::const_iterator itr=value.begin(); itr!=value.end(); itr++){
		std::vector<REAL>::const_iterator itc=(*itr).begin();
		for(int ic=0; ic!=i-1; ic++){
			itc++;		// move to the ith element of itr row
		}
		result.push_back(*itc);
	}
	return result;
}

std::vector<REAL> matrix::getRow(int i){
	if(i<=0 || i>num_row){
		std::cout << "Error: index exceeds!"	<< std::endl;
		exit(-1);
	}
//	std::vector<std::vector<REAL> >::const_iterator itr=value.begin();
	std::vector<std::vector<REAL> >::size_type itr = 0;
	for(int ir=0; ir!=i-1; ir++){
		itr++;
	}
	return value[itr];
}

void matrix::appendRow(std::vector<REAL> row){
	if(row.size() != num_col && num_col != 0){	// num_col != 0 in case that the matrix is empty
		std::cout << "Error: the dimesions do not match!" << std::endl;
		exit(-1);
	}	
	value.push_back(row);
	num_row++;
	num_col = row.size();	// this is important to avoid the case that at first I append a 3 elements row
					// while num_col = 0, num_col should be 3
}

void matrix::appendCol(std::vector<REAL> column){
	if(column.size() != num_row && num_row != 0){
		std::cout << "Error: the dimesions do not match!" << std::endl;
		exit(-1);
	}
	std::vector<REAL> temp_row;
	if(num_row == 0){
		for(std::vector<REAL>::const_iterator itc=column.begin(); itc!=column.end(); itc++){
			temp_row.clear();
			temp_row.push_back(*itc);
			this->appendRow(temp_row);
		}
	}
	else {
		int ir = 0;
		for(std::vector<std::vector<REAL> >::iterator itr=value.begin(); itr!=value.end(); itr++){	// if num_row == 0, then 			means this matrix is empty and this for loop will not be getted in which leads to cannot appendCol to empty matrix
			(*itr).push_back(column[ir]);	// put the ir element of col to the back of itr row
			if(ir>=num_row){
				std::cout << "Error: dimension exceeds in appendCol!" << std::endl;
				exit(-1);
			}
			ir++;
		}
		num_col++;
	}
	num_row = column.size();
}

void matrix::clear(){
	value.clear();
	num_col = 0;
	num_row = 0;
}

matrix matrix::getInvs(){	// at present this code can only get the inverse of a two by two matrix and 1x1 matrix
	if(num_row == 2 && num_col ==2){
		matrix result;
		REAL det;
		REAL a, b, c, d;		// [a b; c d]
		a = value.front().front();
		b = value.front().back();
		c = value.back().front();
		d = value.back().back();
		det = a*d-b*c;
//		if(det == 0){
//			std::cout << "Error: matrix is singular!" << std::endl;
//			exit(-1);
//		}
		std::vector<REAL> temp_row;
		temp_row.clear();	// calculate the first row of result
		temp_row.push_back(d/det);
		temp_row.push_back(-b/det);
		result.appendRow(temp_row);
		temp_row.clear();	// calculate the second row of result
		temp_row.push_back(-c/det);
		temp_row.push_back(a/det);
		result.appendRow(temp_row);
		// return
		result.num_row = num_row;
		result.num_col = num_col;	// get number of rows for the new matrix
		return result;	
	}
	else if(num_row == 1 && num_col == 1){
		matrix result;
		REAL det;
		det = (*this)(1,1);
		det = 1/det;
		std::vector<REAL> temp_row;
		temp_row.push_back(det);
		result.appendRow(temp_row);
		result.num_row = 1;
		result.num_col = 1;
		return result;
	}
	else if(num_row == 3 && num_col ==3){
		matrix result(3,3);
		REAL det;
		REAL a,b,c,d,e,f,g,h,k;
		a = (*this)(1,1);
		b = (*this)(1,2);
		c = (*this)(1,3);
		d = (*this)(2,1);
		e = (*this)(2,2);
		f = (*this)(2,3);
		g = (*this)(3,1);
		h = (*this)(3,2);
		k = (*this)(3,3);
		
		det = a*(e*k-f*h)-b*(k*d-f*g)+c*(d*h-e*g);
		REAL A,B,C,D,E,F,G,H,K;
		A = e*k-f*h;
		B = -(f*g-d*k);
		C = d*h-e*g;
		D = -(b*k-c*h);
		E = a*k-c*g;
		F = -(a*h-b*g);
		G = b*f-c*e;
		H = -(a*f-c*d);
		K = a*e-b*d;

		result(1,1) = A/det;
		result(2,1) = B/det;
		result(3,1) = C/det;
		result(1,2) = D/det;
		result(2,2) = E/det;
		result(3,2) = F/det;
		result(1,3) = G/det;
		result(2,3) = H/det;
		result(3,3) = K/det;

		return result;
	}
	else {
		std::cout << "Sorry: at present this code can only get the inverse of a 3 by 3 matrix!" << std::endl;
		exit(-1);
	}
}

matrix matrix::getTrans(){
	matrix result;
	for(std::vector<std::vector<REAL> >::const_iterator itr=value.begin(); itr!=value.end(); itr++){
		result.appendCol(*itr);
	}
	result.num_row = num_col;
	result.num_col = num_row;	// get number of rows for the new matrix
	return result;
}

REAL matrix::getNorm(){
	if(num_row==1){ 	// for column vector 
		REAL norm = 0;
		for(int ic=1; ic!=num_col+1; ic++){
			norm = norm+((*this)(1,ic))*((*this)(1,ic));
		}
		norm = sqrt(norm);
		return norm;
	}
	else if(num_col==1){	// for row vector
		REAL norm = 0;
		for(int ir=1; ir!=num_row+1; ir++){
			norm = norm+((*this)(ir,1))*((*this)(ir,1));
		}
		norm = sqrt(norm);
		return norm;
	}
	else {
		std::cout << "Sorry: at present we can only get norm of vectors!" << std::endl;
		exit(-1);
	}
}

matrix& matrix::operator = (matrix A){
	// need to delete matrix *this first
	this->clear();
	this->num_row = 0;
	this->num_col = 0;
	for(int ir=0; ir!=A.num_row; ir++){
		this->appendRow(A.getRow(ir+1));
	}
	// get number of rows for the new matrix
	this->num_row = A.num_row;
	this->num_col = A.num_col;
	return *this;
}

REAL& matrix::operator () (int i, int j){
	if(i>num_row || i<1 || j>num_col || j<1){
		std::cout << "Error: index exceeds when ()!" << std::endl;
		std::cout << "i,j in ():\n " << i << ", " << j << std::endl;
		std::cout << "num_row, num_col:\n " << num_row << " " << num_col << std::endl;
		exit(-1);
	}
	std::vector<std::vector<REAL> >::size_type itr = 0;
	for(int ir=0; ir!=i-1; ir++)
		itr++;
	std::vector<REAL>::size_type itc = 0;
	for(int ic=0; ic!=j-1; ic++)
		itc++;
	return value[itr][itc];
}

std::string matrix::print(){
//	std::cout << "Matrix: " << std::endl;
	std::stringstream ss;
	for(std::vector<std::vector<REAL> >::const_iterator itr=value.begin(); itr!=value.end(); itr++){
		for(std::vector<REAL>::const_iterator itc=(*itr).begin(); itc!=(*itr).end(); itc++){
			ss << *itc << " ";
		}
		ss << std::endl;
	}
	return ss.str();
}

//void operator += (matrix A){
//	(*this) = (*this)+A;
//}


// non member functions
matrix operator + (matrix A, matrix B){
	if(A.num_row != B.num_row || A.num_col != B.num_col){
		std::cout << "Error: dimensions do not match!" << std::endl;
		exit(-1);
	}
	matrix result;
	std::vector<REAL> temp_row;	// used to store each row of result matrix
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();	// initialize
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]+B.getRow(ir+1)[ic]);
		}
		result.appendRow(temp_row);
	}
	result.num_row = A.num_row;
	result.num_col = A.num_col;	// get number of rows for the new matrix
	return result;
}

matrix operator + (REAL k, matrix A){
	matrix result;
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]+k);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator + (matrix A, REAL k){
	matrix result;
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]+k);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator - (matrix A, matrix B){
	if(A.num_row != B.num_row || A.num_col != B.num_col){
		std::cout << "Error: dimensions do not match!" << std::endl;
		exit(-1);
	}
	matrix result;	
	std::vector<REAL> temp_row;	// used to store each row of result matrix
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();	// initialize
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]-B.getRow(ir+1)[ic]);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator - (REAL k, matrix A){
	matrix result;
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(k-A.getRow(ir+1)[ic]);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator - (matrix A, REAL k){
	matrix result;
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]-k);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator * (matrix A, matrix B){
	if(A.num_col != B.num_row){
		std::cout << "Error: dimensions do not match!" << std::endl;
		exit(-1);
	}
	matrix result;
	
	REAL temp;
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();	// used to store the ir row of result matrix
		for(int ic=0; ic!=B.num_col; ic++){
			temp = 0;	// initialize, used to store the k element of temp_row
			for(std::vector<REAL>::size_type k=0; k!=A.getRow(ir+1).size(); k++){
				temp += A.getRow(ir+1)[k]*B.getCol(ic+1)[k];
			}
			temp_row.push_back(temp);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = B.num_col;
	return result;
}

matrix operator * (matrix A, REAL k){
	matrix result;
	
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]*k);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator * (REAL k, matrix A){
	matrix result;
	
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]*k);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix operator / (matrix A, REAL k){
	matrix result;
	
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(A.getRow(ir+1)[ic]/k);
		}
		result.appendRow(temp_row);
	}
	// get number of rows for the new matrix
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

matrix expm(matrix A){
	matrix result;
	std::vector<REAL> temp_row;
	for(int ir=0; ir!=A.num_row; ir++){
		temp_row.clear();
		for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
			temp_row.push_back(exp(A.getRow(ir+1)[ic]));
		}
		result.appendRow(temp_row);	
	}
	result.num_row = A.num_row;
	result.num_col = A.num_col;
	return result;
}

REAL trace(matrix A){
	if(A.num_row != A.num_col) {
		std::cout<< "Matrix must be square to calculate trace!" << std::endl;
		exit(-1);
	}
	REAL tra = 0;	// trace
	for(int i=0; i!=A.num_row; i++){
		tra += A(i+1,i+1);
	}
	return tra;
}




}// end of dem


/*
// used to test
using namespace dem;

int main(){
	matrix A;
	std::vector<REAL> temp_row;
	temp_row.push_back(3);
	temp_row.push_back(-4);
	A.appendRow(temp_row);
	std::cout << "should be 3 -4" <<std::endl;
	std::cout << "dimesions of A: " << A.num_row << " " << A.num_col << std::endl;
	std::cout << A.print();
	temp_row.clear();
	temp_row.push_back(-4);
	temp_row.push_back(2);
	A.appendRow(temp_row);
	std::cout << "should be 3 -4; -4 2" << std::endl;
	std::cout << A.print();
	temp_row.clear();
	temp_row.push_back(0);
	temp_row.push_back(1);
	A.appendRow(temp_row);
	std::cout << "should be 3 -4; -4 2; 0 1" << std::endl;
	std::cout << "dimesions of A: " << A.num_row << " " << A.num_col << std::endl;
	std::cout << A.print();
	temp_row.clear();

	std::vector<REAL> temp_col;
	temp_col.clear();
	temp_col.push_back(0);
	temp_col.push_back(1);
	temp_col.push_back(1);
	A.appendCol(temp_col);
	std::cout << "should be 3 -4 0; -4 2 1; 0 1 1" << std::endl;
	std::cout << "dimesions of A: " << A.num_row << " " << A.num_col << std::endl;
	std::cout << A.print();
	
	// test +
	std::cout << "test + begin: " <<std::endl;
	matrix Aplus;
	std::cout << "dimesions of A: " << A.num_row << " " << A.num_col << std::endl;
	std::cout << "dimesions of B: " << A.num_row << " " << A.num_col << std::endl;
	Aplus = A + A;
	std::cout << "should be 6 -8 0; -8 4 2; 0 2 2" << std::endl;
	Aplus.print();
	// test -
	std::cout << "test - begin: " <<std::endl;
	matrix Aplusminus;
	Aplusminus = Aplus - A;
	std::cout << "should be 3 -4 0; -4 2 1; 0 1 1" << std::endl;
	Aplusminus.print();
	// test *
	std::cout << "test * begin: " <<std::endl;
	matrix AA;
	AA = A*A;
	std::cout << "A*A: " << std::endl;
	AA.print();
	matrix AAA;
	AAA = AA*A;
	std::cout << "(A*A)*A: " << std::endl;
	std::cout << AAA.print();

	AAA = A*A*A;
	std::cout << "A*A*A: " << std::endl;
	std::cout << AAA.print();
	// test invs
	std::cout << "test invs begin: " <<std::endl;

	matrix B;
	temp_row.clear();
	temp_row.push_back(8);
	temp_row.push_back(-3);
	B.appendRow(temp_row);

	temp_row.clear();
	temp_row.push_back(-4);
	temp_row.push_back(2);
	B.appendRow(temp_row);
	std::cout << "should be 8 -3; -4 2: " << std::endl;
	std::cout << B.print();
	matrix Binv;
	
	
	Binv = B.getInvs();
	std::cout << "Binvs: " <<std::endl;
	std::cout << Binv.print();
	std::cout << "dimesions of Binv: " << Binv.num_row << " " << Binv.num_col << std::endl;
	matrix BBinv;
	BBinv = Binv*B;
	std::cout << "BBinvs: " << std::endl;
	std::cout << BBinv.print();
	std::cout << "dimesions of BBinv: " << BBinv.num_row << " " << BBinv.num_col << std::endl;
	// test transpose
	std::cout << "test transpose begin: " <<std::endl;
	temp_col.clear();
	temp_col.push_back(0);
	temp_col.push_back(1);
	temp_col.push_back(1);
	A.appendCol(temp_col);
	std::cout << "A: " << std::endl;
	std::cout << A.print();
	std::cout << "dimesions of A: " << A.num_row << " " << A.num_col << std::endl;
	matrix Atran;
	Atran = A.getTrans();
	std::cout << "Atran: " << std::endl;
	std::cout << Atran.print();
	std::cout << "dimesions of Atran: " << Atran.num_row << " " << Atran.num_col << std::endl;
	// test () to see if the elements can be modified
std::cout << "point 1!" << std::endl;
	matrix G(4,1);
//	temp_row.clear();
//	temp_row.push_back(1);
//	G.appendRow(temp_row);
//	G.appendRow(temp_row);
//	G.appendRow(temp_row);
//	G.appendRow(temp_row);

	std::cout << "G: " << std::endl;
	std::cout << G.print();
	std::cout << "dimesions of G: " << G.num_row << " " << G.num_col << std::endl;
	G(2,1) = 1;
	std::cout << "G: " << std::endl;
	std::cout << G.print();
	std::cout << "dimesions of G: " << G.num_row << " " << G.num_col << std::endl;
	std::cout << "2*G: " << std::endl;
	std::cout << (2*G).print();
	// test getNorm()
	std::cout << "G norm: " << G.getNorm() << std::endl;
	// test exp
	std::cout << "2G exp: " <<  std::endl;
	std::cout << expm(G).print();
	// tesp -/+
	std::cout << "G+1: " << std::endl;
	std::cout << (G+1).print();
	std::cout << "1+G: " << std::endl;
	std::cout << (1+G).print();
	std::cout << "G-1: " << std::endl;
	std::cout << (G-1).print();
	std::cout << "1-G: " << std::endl;
	G -= 1;
	G -= 2*G;
	std::cout << (G).print();
	// test matrix(3,4)
	matrix F(3,4);
	std::cout << "F: " << std::endl;
	std::cout << F.print();
	std::cout << "dimension: " << F.num_row << " " << F.num_col << std::endl;
	F(2,4) = 24;
	std::cout << "F: " << std::endl;
	std::cout << (2+F).print();
	std::cout << "dimension: " << F.num_row << " " << F.num_col << std::endl;
	// test matrix += 4
//	2+F;
	std::cout << "F: " << std::endl;
//	std::cout << (2.0+F).print();
//	std::cout << "dimension: " << F.num_row << " " << F.num_col << std::endl;
	F *= G;
//	F += (2*G);
	std::cout << "F: " << std::endl;
	std::cout << F.print();
	std::cout << "dimension: " << F.num_row << " " << F.num_col << std::endl;
	// test /
	F = F/2;
	std::cout << "F: " << std::endl;
	std::cout << F.print();
	std::cout << "dimension: " << F.num_row << " " << F.num_col << std::endl;

}
*/


