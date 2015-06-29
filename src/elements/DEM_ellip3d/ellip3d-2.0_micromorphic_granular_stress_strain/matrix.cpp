#include "matrix.h"
#include <cstdlib>
#include <cmath>
// written on Feb 13, 2013
namespace dem {

matrix::matrix(int i, int j){
	if(i<0 || j<0){
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
//		calcDimensions();

	}
}

std::vector<REAL> matrix::getCol(int i) const{
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

std::vector<REAL> matrix::getRow(int i) const{
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


void matrix::calcDimensions(){
	
    std::vector<std::vector<REAL> >::size_type row = value.size();
    std::vector<REAL>::size_type col = value[0].size();

    num_row = row;
    num_col = col;

} // calcDimensions()


void matrix::appendRow(const std::vector<REAL> & row){
	if(row.size() != num_col && num_col != 0){	// num_col != 0 in case that the matrix is empty
		std::cout << "Error: the dimesions do not match!" << std::endl;
		exit(-1);
	}	
	value.push_back(row);
	num_row++;
	num_col = row.size();	// this is important to avoid the case that at first I append a 3 elements row

	calcDimensions();
					// while num_col = 0, num_col should be 3
}

void matrix::appendCol(const std::vector<REAL> & column){
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
	calcDimensions();
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
		result.calcDimensions();
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
		result.calcDimensions();
		return result;
	}
	else if(num_row == 3 && num_col == 3){
		double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13;
		double z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25;	

		z1 = (*this)(1,1);
		z2 = (*this)(1,2);
		z3 = (*this)(1,3);
		z4 = (*this)(2,1);
		z5 = (*this)(2,2);
		z6 = (*this)(2,3);
		z7 = (*this)(3,1);
		z8 = (*this)(3,2);
		z9 = (*this)(3,3);
		z10 = -z2*z4;
		z11 = z3*z4;
		z12 = z1*z5;
		z13 = -z3*z5;
		z14 = -z1*z6;
		z15 = z2*z6;
		z16 = z13*z7;
		z17 = z15*z7;
		z18 = z2*z7;
		z19 = -z3*z7;
		z20 = -z5*z7;
		z7 = z6*z7;
		z21 = -z1*z8;
		z22 = z11*z8;
		z23 = z14*z8;
		z3 = z3*z8;
		z24 = z4*z8;
		z6 = -z6*z8;
		z1 = z1*z9;
		z8 = z10*z9;
		z25 = z12*z9;
		z2 = -z2*z9;
		z4 = -z4*z9;
		z5 = z5*z9;
		z9 = z10 + z12;
		z10 = z11 + z14;
		z11 = z13 + z15;
		z12 = z18 + z21;
		z13 = z20 + z24;
		z1 = z1 + z19;
		z8 = z16 + z17 + z22 + z23 + z25 + z8;
		z2 = z2 + z3;
		z3 = z4 + z7;
		z4 = z5 + z6;
		z5 = 1.0/z8;
		z6 = z5*z9;
		z7 = z10*z5;
		z8 = z11*z5;
		z9 = z12*z5;
		z10 = z13*z5;
		z1 = z1*z5;
		z2 = z2*z5;
		z3 = z3*z5;
		z4 = z4*z5;

		//{{z4, z2, z8},
		// {z3, z1, z7},
		// {z10, z9, z6}}
		
		matrix result(3,3);
		result(1,1) = z4;
		result(1,2) = z2;
		result(1,3) = z8;
		result(2,1) = z3;
		result(2,2) = z1;
		result(2,3) = z7;
		result(3,1) = z10;
		result(3,2) = z9;
		result(3,3) = z6;
		result.calcDimensions();
		return result;
	}
	else {
		std::cout << "Sorry: at present this code can only get the inverse of a 2 by 2 matrix and 3x3 matrix!" << std::endl;
		exit(-1);
	}

}

matrix matrix::getTrans() const{
	matrix result;
	for(std::vector<std::vector<REAL> >::const_iterator itr=value.begin(); itr!=value.end(); itr++){
		result.appendCol(*itr);
	}
	result.num_row = num_col;
	result.num_col = num_row;	// get number of rows for the new matrix
	result.calcDimensions();
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

void matrix::LU(matrix & L, matrix & U) const{

    if(num_col!=num_row){
	std::cout << "should be square matrix in LU decomposion..." << std::endl;
	exit(-1);
    }

    matrix a = (*this);

    L = zeros(num_row, num_col);
    U = zeros(num_row, num_col);

//#pragma omp parallel shared(a, L, U)
//{
//    #pragma omp for schedule(static, 10)
    for(int k=1; k<num_row+1; ++k){
	L(k,k) = 1;
	for(int i=k+1; i<num_row+1; ++i){
	    L(i, k) = a(i, k)/a(k,k);
	    // a(i,k) = L(i,k);
	    for(int j=k+1; j<num_row+1; ++j){
		a(i,j) = a(i,j)-L(i,k)*a(k,j);
	    }
//	    #pragma omp flush(a)
	}
//	#pragma omp nowait
	for(int j=k; j<num_row+1; ++j){
	    U(k,j) = a(k,j);
	}
    }

//} // end parallel

} // LU()


REAL matrix::getVal(int i, int j) const {

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

} // getVal()


matrix& matrix::operator = (const matrix &A){
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
	calcDimensions();
	return *this;
}

REAL& matrix::operator () (int i, int j) {
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

std::string matrix::print() const{
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
matrix operator + (const matrix & A, const matrix & B){
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
	result.calcDimensions();
	return result;
}

matrix operator + (REAL k, const matrix & A){
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
	result.calcDimensions();
	return result;
}

matrix operator + (const matrix & A, REAL k){
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
	result.calcDimensions();
	return result;
}

matrix operator - (const matrix & A, const matrix & B){
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
	result.calcDimensions();
	return result;
}

matrix operator - (REAL k, const matrix & A){
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
	result.calcDimensions();
	return result;
}

matrix operator - (const matrix & A, REAL k){
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
	result.calcDimensions();
	return result;
}

matrix operator * (const matrix & A, const matrix & B){
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
	result.calcDimensions();
	return result;
}

matrix operator * (const matrix & A, REAL k){
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
	result.calcDimensions();
	return result;
}

matrix operator * (REAL k, const matrix & A){
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
	result.calcDimensions();
	return result;
}

matrix operator / (const matrix & A, REAL k){
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
	result.calcDimensions();
	return result;
}


matrix operator % (matrix & A, matrix & r){ // left division, "\"

    matrix x;
    matrixEqnSolver(x, A, r);

    return x;

} // % 


matrix expm(const matrix & A){
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
	result.calcDimensions();
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


bool isnan_mat(matrix &mat){
    
    bool is = false;
    for(int i=1; i<mat.num_row+1; ++i){
	for(int j=1; j<mat.num_col+1; ++j){
	    if( std::isnan(mat(i,j)) )
		is = true;
	}
    }

    return is;

} // isnan_mat()


int size(const matrix & mat, int n){
    if(n!=1 && n!=2){
	std::cout << "Only can get number of rows and number of columns with size()." << std::endl;
	exit(-1);
    }
    else if(n==1){
	return mat.num_row;
    }
    else{
	return mat.num_col;
    }
} // size()


matrix ones(int i, int j){

    matrix mat;

    if(i<0 || j<0){
	std::cout << "Error: the dimensions of matrix should be positive!" << std::endl;
	exit(-1);	
    }
    else if(i==0 || j==0){
	return mat;
    }
    else{
	mat.num_row = 0;
	mat.num_col = 0;	// to avoid errors in appendRow()
	std::vector<REAL> temp_vec;
	for(int ir=0; ir!=i; ir++){
	    temp_vec.clear();	// begin a new row			
 	    for(int ic=0; ic!=j; ic++){
		temp_vec.push_back(1);
	    }
	    mat.appendRow(temp_vec);
	}
	mat.num_row = i;
	mat.num_col = j;
//	mat.calcDimensions();
    }

    return mat; 

} // ones()


matrix zeros(int i, int j){

    matrix mat(i,j);
    return mat;

} // zeros()


REAL max(matrix & mat){

    REAL max = mat(1,1);
    for(std::vector<std::vector<REAL> >::const_iterator itr=mat.value.begin(); itr!=mat.value.end(); itr++){
	for(std::vector<REAL>::const_iterator itc=(*itr).begin(); itc!=(*itr).end(); itc++){
	    if(max<(*itc))
		max = (*itc);
	}
    }

    return max;

} // max()


REAL min(matrix & mat){

    REAL min = mat(1,1);
    for(std::vector<std::vector<REAL> >::const_iterator itr=mat.value.begin(); itr!=mat.value.end(); itr++){
	for(std::vector<REAL>::const_iterator itc=(*itr).begin(); itc!=(*itr).end(); itc++){
	    if(min>(*itc))
		min = (*itc);
	}
    }

    return min;

} // min()


matrix abs(matrix & A){
	
    matrix result;
    std::vector<REAL> temp_row;
    for(int ir=0; ir!=A.num_row; ir++){
	temp_row.clear();
	for(std::vector<REAL>::size_type ic=0; ic!=A.getRow(ir+1).size(); ic++){
  	    temp_row.push_back(fabs(A.getRow(ir+1)[ic]));
	}
	result.appendRow(temp_row);	
    }
    result.num_row = A.num_row;
    result.num_col = A.num_col;
    result.calcDimensions();
    return result;

} // abs()


int length(const matrix & mat){

    int num;
    if(mat.num_row==0 || mat.num_col==0){
	num = 0;
    }
    else {
    	num = mat.num_row;
    	if(num < mat.num_col)
	    num = mat.num_col;
    }
    return num;

} // length()


REAL norm(matrix & mat){

    REAL val = mat.getNorm();

    return val;

} // norm()


REAL det(matrix &mat){

    if(mat.num_row==2 && mat.num_col==2){ // (2 x 2)

	return mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1);
    }
    else if(mat.num_row==1 && mat.num_col==1){ // (1 x 1)

	return mat(1,1);
    }
    else if(mat.num_row==3 && mat.num_col==3){ // (3 x 3)

	return  mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))
	      - mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1))
	      + mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1));

    }
    else{
	std::cout << "Matrix dimensions problem in det()..." << std::endl;
	exit(-1);
    }

} // det()


matrix linspace(REAL start, REAL stop, int num){

    matrix result(1,num);

    if(num<=0){
	std::cout << "Number should be positive in linspace()..." << std::endl;
	exit(-1);
    }	
    else if(num==1){
	result(1,1) = stop;
    }
    else{

	REAL step = (stop-start)/double((num-1));
	for(int i=1; i<num; ++i){
	    result(1,i) = start+step*(i-1);
	}
	result(1,num) = stop;
    }
    result.calcDimensions();
    return result;


} // linspace()


void matrixEqnSolver(matrix &x, matrix&A, matrix&r){

    if(A.num_row!=r.num_row || r.num_col!=1){
	std::cout << "Dimensions do not match in matrixEqnSolver()!" << std::endl;
	exit(-1); 
    }

    matrix L, U;
    A.LU(L, U);

    matrix y(A.num_row, 1);
    for(int i=1; i<A.num_row+1; ++i){
	REAL left = 0;
	for(int j=1; j<i; ++j){
	    left += L(i,j)*y(j,1);
	}
	y(i,1) = (r(i,1)-left)/L(i,i);
    }

    x = zeros(r.num_row,1);

    for(int i=A.num_row; i>0; --i){
	REAL left = 0;
	for(int j=i+1; j<A.num_col+1; ++j){
	    left += U(i,j)*x(j,1);
	}
	x(i,1) = (y(i,1)-left)/U(i,i);
    }

} // matrixEqnSolver()


}// end of memFluid


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


