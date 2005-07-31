/* $Id: dSymMatrixT.cpp,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (03/03/1997)                                          */

#include "dSymMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "Constants.h"
#include "dMatrixT.h"

const double Pi = acos(-1);
inline double d_sign(double a, double b)
{
	double x;
	x = (a >= 0 ? a : -a);
	return (b >= 0 ? x : -x);
};

/* constructor */
dSymMatrixT::dSymMatrixT(void): fNumSD(0) { }
dSymMatrixT::dSymMatrixT(int nsd) { Allocate(nsd); }
dSymMatrixT::dSymMatrixT(int nsd, double* array) { Set(nsd,array); }
dSymMatrixT::dSymMatrixT(const dSymMatrixT& source): fNumSD(0)
{
	operator=(source);
}

/* allocate an reduced matrix for the given spatial dimension */
void dSymMatrixT::Allocate(int nsd)
{
	/* check  */
	fNumSD = nsd;
	if (fNumSD < 1 && fNumSD > 3) throw eGeneralFail;

	/* inherited */
	dArrayT::Allocate(NumValues(fNumSD));
}

/* set fields */
void dSymMatrixT::Set(int nsd, double* array)
{
	fNumSD = nsd;
	if (fNumSD < 1 && fNumSD > 3) throw eGeneralFail;
	
	/* inherited */
	dArrayT::Set(NumValues(fNumSD), array);
}

//DEV - this is bad
#if 0
void dSymMatrixT::SetNSD(int nsd)
{
	fNumSD = nsd;
	if (fNumSD < 1 && fNumSD > 3) throw eGeneralFail;
}
#endif

/* accessor */
double& dSymMatrixT::operator()(int row, int col) const
{
/* rigorous range checking */
#if __option (extended_errorcheck)
	if (row < 0 || row >= fNumSD ||
	    col < 0 || col >= fNumSD) throw eOutOfRange;
#endif

	int map2D[2][2] = {{0,2},
	                   {2,1}};

	int map3D[3][3] = {{0,5,4},
	                   {5,1,3},
	                   {4,3,2}};

	if (fNumSD == 0) throw eGeneralFail; //not configured
	if (fNumSD == 2)
		return fArray[map2D[row][col]];		
	else
		return fArray[map3D[row][col]];		
}

void dSymMatrixT::ExpandIndex(int nsd, int dex, int& dex_1, int& dex_2)
{
#if __option(extended_errorcheck)
	/* consistency check */
	if (dex >= NumValues(nsd)) throw eOutOfRange;
#endif	

	int  map_2D[6] = {0,0,1,1,0,1};
	int map_3D[18] = {0,0,1,1,2,2,
	                  1,2,0,2,0,1};
	int* map = (nsd == 2) ? map_2D : map_3D;

	int* p = map + 2*dex;
	dex_1 = p[0];
	dex_2 = p[1];
}
	
/* I/O operators */
ostream& operator<<(ostream& out, const dSymMatrixT& array)
{
	if (array.fNumSD == 0) throw eGeneralFail;

	int d_width = OutputWidth(out, array.fArray);

	if (array.fNumSD == 2)
	{
		out << setw(d_width) << array.fArray[0];
		out << setw(d_width) << array.fArray[2] << '\n';
		out << setw(d_width) << array.fArray[2];
		out << setw(d_width) << array.fArray[1] << '\n';
	}
	else if (array.fNumSD == 3)
	{
		out << setw(d_width) << array.fArray[0];
		out << setw(d_width) << array.fArray[5];
		out << setw(d_width) << array.fArray[4] << '\n';
		out << setw(d_width) << array.fArray[5];
		out << setw(d_width) << array.fArray[1];
		out << setw(d_width) << array.fArray[3] << '\n';
		out << setw(d_width) << array.fArray[4];
		out << setw(d_width) << array.fArray[3];
		out << setw(d_width) << array.fArray[2] << '\n';
	}
	else
		out << setw(d_width) << array.fArray[0] << '\n';

	return out;
}

istream& operator>>(istream& in, dSymMatrixT& array)
{
	if (array.fNumSD == 2)
	{
		double a21;

		in >> array.fArray[0] >> array.fArray[2];
		in >> a21             >> array.fArray[1];
		
		if (a21 != array.fArray[2]) throw eBadInputValue;
	}
	else if (array.fNumSD == 3)
	{
		double a21, a31, a32;

		in >> array.fArray[0] >> array.fArray[5] >> array.fArray[4];
		in >> a21             >> array.fArray[1] >> array.fArray[3];
		in >> a31             >> a32             >> array.fArray[2];
		
		if (a21 != array.fArray[5] ||
		    a31 != array.fArray[4] ||
		    a32 != array.fArray[3]) throw eBadInputValue;
	}
	else if (array.fNumSD == 1)
		in >> array.fArray[0];
	else
		throw eGeneralFail;

	return in;
}


/* return eigenvalues and eigenvectors (in columns) */
void dSymMatrixT::PrincipalValues(dArrayT& val) const // will get phased out
{
	/* redirect */
	Eigenvalues(val, true);
}

void dSymMatrixT::Eigenvalues(dArrayT& val, bool sort_descending) const
{
	if (fNumSD == 1)
		val[0] = fArray[0];
	else if (fNumSD == 2)
	{
	double mean = 0.5*(fArray[0] + fArray[1]);
	double diff = 0.5*(fArray[0] - fArray[1]);
	double r    = sqrt(fArray[2]*fArray[2] + diff*diff);

	val[0] = mean + r;
	val[1] = mean - r;
	}
	else if (fNumSD == 3)
	{
		Eigenvalues3D_Cardano(val);

		//Eigenvalues3D(val, sort_descending);

//TEMP - iterative eigs not guaranteed to be in descending order
//       but it's not clear what assumes the eigs will be ordered
#if 0
		dArrayT tmp(3);		
		Eigenvalues3D(tmp);
		for (int i = 0; i < 3; i++)
			if (fabs(val[0] - tmp[0]) > kSmall) throw eGeneralFail;
#endif
	}
	else
		throw eGeneralFail;
}

void dSymMatrixT::Eigensystem(dArrayT& val, dMatrixT& vec, bool sort_descending) const
{
	if (fNumSD == 1)
	{
		val[0] = fArray[0];
		vec[0] = 1.0;
	}
	else if (fNumSD == 2)
	{
	double mean = 0.5*(fArray[0] + fArray[1]);
	double diff = 0.5*(fArray[0] - fArray[1]);
	double r    = sqrt(fArray[2]*fArray[2] + diff*diff);

		/* eigenvalues (sorted) */
	val[0] = mean + r;
	val[1] = mean - r;
	
	/* eigenvectors */
	double t;
	if (fabs(val[0] - val[1]) < kSmall)
		t = 0.0;
	else
	{
		double arg = diff/r;
			t = (fArray[2] > 0) ?
			0.5*acos(arg) :
			0.5*(acos(-arg) + Pi); // angle to largest eigenvalue
	}

		vec[0] = cos(t);
		vec[1] = sin(t);
		
		vec[2] =-vec[1];
		vec[3] = vec[0];
	}
	else if (fNumSD == 3)
		Eigensystem3D(val, vec, sort_descending);
	else
		throw eGeneralFail;
}

/*
* Returns the scalar inner product of a tensor with itself
* define as:
*
*             T:T = T_ij T_ij
*
*/
double dSymMatrixT::ScalarProduct(void) const
{
#if __option(extended_errorcheck)
	if (fNumSD == 0) throw eGeneralFail;
#endif	

	//return( 2.0*Invariant2() + pow(Trace(),2) );
	if (fNumSD == 2)
		return fArray[0]*fArray[0] +
		       fArray[1]*fArray[1] +
		   2.0*fArray[2]*fArray[2];
	else if (fNumSD == 3)
		return fArray[0]*fArray[0] +
		       fArray[1]*fArray[1] +
		       fArray[2]*fArray[2] +
		  2.0*(fArray[3]*fArray[3] +
		       fArray[4]*fArray[4] +
		       fArray[5]*fArray[5]);
	else
		return fArray[0]*fArray[0];
}

/* returns the magnitude of the second invariant of the reduced index
* tensor.  For 2D, the out-of-plane components are assumed to be zero,
* ie if *this represents stress, then plane stress is assumed and if
* *this is strain then plane strain is assumed. */
double dSymMatrixT::Invariant2(void) const
{
#if __option(extended_errorcheck)
	if (fNumSD == 0) throw eGeneralFail;
#endif	

	if (fNumSD == 2)
		return fArray[2]*fArray[2] -
		       fArray[0]*fArray[1];
	else if (fNumSD == 3)
		return fArray[3]*fArray[3] +
			   fArray[4]*fArray[4] +
			   fArray[5]*fArray[5] -
			   fArray[0]*fArray[1] -
			   fArray[1]*fArray[2] -
			   fArray[2]*fArray[0];
	else
		return 0.0;
}

/* matrix determinant */
double dSymMatrixT::Det(void) const
{
#if __option(extended_errorcheck)
	if (fNumSD == 0) throw eGeneralFail;
#endif	

	if (fNumSD == 2)
		return fArray[0]*fArray[1] - fArray[2]*fArray[2];
	else if (fNumSD == 3)
		return fArray[0]*(fArray[1]*fArray[2] - fArray[3]*fArray[3])
			 - fArray[5]*(fArray[5]*fArray[2] - fArray[3]*fArray[4])
			 + fArray[4]*(fArray[5]*fArray[3] - fArray[1]*fArray[4]);
	else
		return fArray[0];
}

/* returns the trace of *this.  For 2D, the out-of-plane
* component is assumed to be zero, ie if *this represents stress,
* then plane stress is assumed and if *this is strain then plane
* strain is assumed. */
double dSymMatrixT::Trace(void) const
{
#if __option(extended_errorcheck)
	if (fNumSD == 0) throw eGeneralFail;
#endif	

	if (fNumSD == 2)
		return fArray[0] + fArray[1];
	else if (fNumSD == 3)
		return fArray[0] + fArray[1] + fArray[2];
	else
		return fArray[0];
	
}

/* identity operations */
void dSymMatrixT::PlusIdentity(double value)
{
#if __option(extended_errorcheck)
	if (fNumSD == 0) throw eGeneralFail;
#endif	

	if (fNumSD == 2)
	{
		fArray[0] += value;
		fArray[1] += value;
	}
	else if (fNumSD == 3)
	{
		fArray[0] += value;
		fArray[1] += value;
		fArray[2] += value;
	}
	else
		fArray[0] += value;
}

dSymMatrixT& dSymMatrixT::Identity(double value)
{
#if __option(extended_errorcheck)
	if (fNumSD == 0) throw eGeneralFail;
#endif	

	if (fNumSD == 2)
	{
		fArray[0] = value;
		fArray[1] = value;
		fArray[2] = 0.0;
	}
	else if (fNumSD == 3)
	{
		fArray[0] = value;
		fArray[1] = value;
		fArray[2] = value;
		fArray[3] = 0.0;
		fArray[4] = 0.0;
		fArray[5] = 0.0;		
	}
	else
		fArray[0] = value;
		
	return *this;
}

/* returns the deviatoric part of *this.
*
* Note: This function should only be called when working in
* 3D, ie. a vector length 6, since there's no way to enforce
* Trace(Deviatoric{vector)) = 0, in 2D. */
dSymMatrixT& dSymMatrixT::Deviatoric(const dSymMatrixT& tensor)
{
	/* dimension checks */
	if (fNumSD != 3) throw eGeneralFail;

	double* pLHS = fArray;
	double* pRHS = tensor.fArray;
	double trace = (pRHS[0] + pRHS[1] + pRHS[2])/3.0;

	*pLHS++ = *pRHS++ - trace;
	*pLHS++ = *pRHS++ - trace;
	*pLHS++ = *pRHS++ - trace;
	*pLHS++ = *pRHS++;
	*pLHS++ = *pRHS++;
	*pLHS   = *pRHS;

	return *this;
}

/* matrix inverse */
dSymMatrixT& dSymMatrixT::Inverse(const dSymMatrixT& matrix)
{
	if (fNumSD == 2)
	{
		double a0 = matrix.fArray[0];
		double a1 = matrix.fArray[1];
		double a2 = matrix.fArray[2];
		double det = a0*a1 - a2*a2;

		fArray[0] = a1/det;
		fArray[1] = a0/det;
		fArray[2] =-a2/det;	
	}
	else if (fNumSD == 3)
	{
		double a0 = matrix.fArray[0];
		double a1 = matrix.fArray[1];
		double a2 = matrix.fArray[2];
		double a3 = matrix.fArray[3];
		double a4 = matrix.fArray[4];
		double a5 = matrix.fArray[5];

		double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
		double z13, z14, z15;

		z1 = a0*a1;
		z2 = a0*a2;
		z3 = a1*a2;
		z4 = -a0*a3;
		z5 = a3*a3;
		z6 = -a1*a4;
		z7 = a3*a4;
		z8 = a4*a4;
		z9 = -a2*a5;
		z10 = a3*a5;
		z11 = a4*a5;
		z12 = a5*a5;
		z13 = a2*z1;
		z14 = 2.*a5*z7;
		z12 = -z12;
		z15 = a2*z12;
		z4 = z11 + z4;
		z5 = -z5;
		z11 = a0*z5;
		z6 = z10 + z6;
		z8 = -z8;
		z10 = a1*z8;
		z7 = z7 + z9;
		z1 = z1 + z12;
		z2 = z2 + z8;
		z8 = z10 + z11 + z13 + z14 + z15;
		z3 = z3 + z5;
		z5 = 1./z8;
		z4 = z4*z5;
		z6 = z5*z6;
		z7 = z5*z7;
		z1 = z1*z5;
		z2 = z2*z5;
		z3 = z3*z5;	
	
		// {{z3, z7, z6},
		//  {z7, z2, z4},
		//  {z6, z4, z1}}
		fArray[0] = z3;
		fArray[1] = z2;
		fArray[2] = z1;
		fArray[3] = z4;
		fArray[4] = z6;
		fArray[5] = z7;
	}
	else
		fArray[0] = 1.0/matrix.fArray[0];
		
	return *this;
}

/* reduced index <-> matrix transformations */	
void dSymMatrixT::ToMatrix(dMatrixT& matrix) const
{
	/* must be square and 2D or 3D */
#if __option(extended_errorcheck)
	if (matrix.Rows() != matrix.Cols() ||
		matrix.Rows() != fNumSD) throw eSizeMismatch;
#endif

	double* pmat = matrix.Pointer();

	if (fNumSD == 2)
	{
		pmat[0] = fArray[0];			//1,1
		pmat[3] = fArray[1];			//2,2
		pmat[1] = pmat[2] = fArray[2];	//1,2
	}
	else if (fNumSD == 3)
	{
		pmat[0] = fArray[0];			//1,1
		pmat[4] = fArray[1];			//2,2
		pmat[8] = fArray[2];			//3,3
		pmat[7] = pmat[5] = fArray[3];	//2,3
		pmat[6] = pmat[2] = fArray[4];	//1,3
		pmat[3] = pmat[1] = fArray[5];	//1,2
	}
	else if (fNumSD == 1)
		pmat[0] = fArray[0];			//1,1
	else
		throw eGeneralFail;
}

/* take the symmetric part of the matrix. Returns a reference to *this */
dSymMatrixT& dSymMatrixT::Symmetrize(const dMatrixT& matrix)
{
	/* dimension check */
	if (fNumSD != matrix.Rows() ||
	    fNumSD != matrix.Cols()) throw eSizeMismatch;

	double* pmat = matrix.Pointer();

	if (fNumSD == 2)
	{
		fArray[0] = pmat[0];
		fArray[1] = pmat[3];

		fArray[2] = 0.5*(pmat[1] + pmat[2]);	
	}
	else if (fNumSD == 3)
	{
		fArray[0] = pmat[0];
		fArray[1] = pmat[4];
		fArray[2] = pmat[8];

		fArray[3] = 0.5*(pmat[5] + pmat[7]);	
		fArray[4] = 0.5*(pmat[2] + pmat[6]);
		fArray[5] = 0.5*(pmat[1] + pmat[3]);
	}
	else
		fArray[0] = pmat[0];

	return *this;
}

dSymMatrixT& dSymMatrixT::FromMatrix(const dMatrixT& matrix)
{
	/* must be square and 2D or 3D */
#if __option(extended_errorcheck)
	if (matrix.Rows() != matrix.Cols() ||
		matrix.Rows() != fNumSD) throw eSizeMismatch;
#endif

	double* pmat = matrix.Pointer();

	if (fNumSD == 2)
	{
		fArray[0] = pmat[0];	//1,1
		fArray[1] = pmat[3];	//2,2
		fArray[2] = pmat[2];	//1,2
	}
	else if (fNumSD == 3)
	{
		fArray[0] = pmat[0];	//1,1
		fArray[1] = pmat[4];	//2,2
		fArray[2] = pmat[8];	//3,3
		fArray[3] = pmat[7];	//2,3
		fArray[4] = pmat[6];	//1,3
		fArray[5] = pmat[3];	//1,2
	}
	else if (fNumSD == 1)
		fArray[0] = pmat[0];	//1,1
	else
		throw eGeneralFail;
	
	return *this;
}

/* 2D <-> 3D translations.  Return references to *this */
dSymMatrixT& dSymMatrixT::ExpandFrom2D(const dSymMatrixT& vec2D) /* assumed plane strain */
{
/* dimension check */
#if __option (extended_errorcheck)
	if (fNumSD != 3 ||
	    vec2D.fNumSD != 2) throw eSizeMismatch;
#endif

	/* translation */
	fArray[0] = vec2D.fArray[0];
	fArray[1] = vec2D.fArray[1];
	fArray[2] = 0.0;
	fArray[3] = 0.0;
	fArray[4] = 0.0;
	fArray[5] = vec2D.fArray[2];
	
	return *this;
}

dSymMatrixT& dSymMatrixT::ReduceFrom3D(const dSymMatrixT& vec3D)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fNumSD != 2 ||
	    vec3D.fNumSD != 3) throw eSizeMismatch;
#endif

	/* translation */
	fArray[0] = vec3D.fArray[0];
	fArray[1] = vec3D.fArray[1];
	fArray[2] = vec3D.fArray[5];
	
	return *this;
}

/* symmetric rank-4 - rank-2 contraction */
void dSymMatrixT::A_ijkl_B_kl(const dMatrixT& A, const dSymMatrixT& B)
{
#if __option(extended_errorcheck)
	if (A.Cols() != NumValues(B.Rows()) ||
	    A.Rows() != NumValues(Rows())   ||
	      Rows() != B.Rows()) throw eSizeMismatch;
#endif

	double* pA = A.Pointer();
	double* pB = B.Pointer();
	if (fNumSD == 2)
	{
		fArray[0] = pA[0]*pB[0] + pA[3]*pB[1] + 2*pA[6]*pB[2];
		fArray[1] = pA[1]*pB[0] + pA[4]*pB[1] + 2*pA[7]*pB[2];
		fArray[2] = pA[2]*pB[0] + pA[5]*pB[1] + 2*pA[8]*pB[2];
	}
	else if (fNumSD == 3)
	{
		fArray[0] = pA[0]*pB[0] + pA[6]*pB[1] + pA[12]*pB[2] + 2*pA[18]*pB[3] + 2*pA[24]*pB[4] + 2*pA[30]*pB[5];
		fArray[1] = pA[1]*pB[0] + pA[7]*pB[1] + pA[13]*pB[2] + 2*pA[19]*pB[3] + 2*pA[25]*pB[4] + 2*pA[31]*pB[5];
		fArray[2] = pA[2]*pB[0] + pA[8]*pB[1] + pA[14]*pB[2] + 2*pA[20]*pB[3] + 2*pA[26]*pB[4] + 2*pA[32]*pB[5];
		fArray[3] = pA[3]*pB[0] + pA[9]*pB[1] + pA[15]*pB[2] + 2*pA[21]*pB[3] + 2*pA[27]*pB[4] + 2*pA[33]*pB[5];
		fArray[4] = pA[4]*pB[0] + pA[10]*pB[1] + pA[16]*pB[2] + 2*pA[22]*pB[3] + 2*pA[28]*pB[4] + 2*pA[34]*pB[5];
		fArray[5] = pA[5]*pB[0] + pA[11]*pB[1] + pA[17]*pB[2] + 2*pA[23]*pB[3] + 2*pA[29]*pB[4] + 2*pA[35]*pB[5];
	}
	else if (fNumSD == 1)
		fArray[0] = pA[0]*pB[0];
	else
		throw eGeneralFail;
}

void dSymMatrixT::A_ijkl_B_ij(const dMatrixT& A, const dSymMatrixT& B)
{
#if __option(extended_errorcheck)
	if (A.Rows() != NumValues(B.Rows()) ||
	    A.Cols() != NumValues(Rows())   ||
	      Rows() != B.Rows()) throw eSizeMismatch;
#endif

	double* pA = A.Pointer();
	double* pB = B.Pointer();
	if (fNumSD == 2)
	{
		fArray[0] = pA[0]*pB[0] + pA[1]*pB[1] + 2*pA[2]*pB[2];
		fArray[1] = pA[3]*pB[0] + pA[4]*pB[1] + 2*pA[5]*pB[2];
		fArray[2] = pA[6]*pB[0] + pA[7]*pB[1] + 2*pA[8]*pB[2];
	}
	else if (fNumSD == 3)
	{
		fArray[0] = pA[0]*pB[0] + pA[1]*pB[1] + pA[2]*pB[2] + 2*pA[3]*pB[3] + 2*pA[4]*pB[4] + 2*pA[5]*pB[5];
		fArray[1] = pA[6]*pB[0] + pA[7]*pB[1] + pA[8]*pB[2] + 2*pA[9]*pB[3] + 2*pA[10]*pB[4] + 2*pA[11]*pB[5];
		fArray[2] = pA[12]*pB[0] + pA[13]*pB[1] + pA[14]*pB[2] + 2*pA[15]*pB[3] + 2*pA[16]*pB[4] + 2*pA[17]*pB[5];
		fArray[3] = pA[18]*pB[0] + pA[19]*pB[1] + pA[20]*pB[2] + 2*pA[21]*pB[3] + 2*pA[22]*pB[4] + 2*pA[23]*pB[5];
		fArray[4] = pA[24]*pB[0] + pA[25]*pB[1] + pA[26]*pB[2] + 2*pA[27]*pB[3] + 2*pA[28]*pB[4] + 2*pA[29]*pB[5];
		fArray[5] = pA[30]*pB[0] + pA[31]*pB[1] + pA[32]*pB[2] + 2*pA[33]*pB[3] + 2*pA[34]*pB[4] + 2*pA[35]*pB[5];
	}
	else if (fNumSD == 1)
		fArray[0] = pA[0]*pB[0];
	else
		throw eGeneralFail;
}

/* symmetric rank-(2-4-2) contraction with this */
double dSymMatrixT::B_ij_A_ijkl_B_kl(const dMatrixT& A) const
{
#if __option(extended_errorcheck)
	if (A.Rows() != A.Cols() ||
	    A.Rows() != NumValues(Rows())) throw eSizeMismatch;
#endif

	double* pA = A.Pointer();
	double* pB = Pointer();
	if (fNumSD == 2)
		return pB[0]*(pA[0]*pB[0] + pA[1]*pB[1] + 2*pA[2]*pB[2]) +
pB[1]*(pA[3]*pB[0] + pA[4]*pB[1] + 2*pA[5]*pB[2]) +
2*pB[2]*(pA[6]*pB[0] + pA[7]*pB[1] + 2*pA[8]*pB[2]);
	else if (fNumSD == 3)
		return pB[0]*(pA[0]*pB[0] + pA[1]*pB[1] + pA[2]*pB[2] +
		            2*pA[3]*pB[3] + 2*pA[4]*pB[4] + 2*pA[5]*pB[5]) +
pB[1]*(pA[6]*pB[0] + pA[7]*pB[1] + pA[8]*pB[2] +
2*pA[9]*pB[3] + 2*pA[10]*pB[4] + 2*pA[11]*pB[5]) +
pB[2]*(pA[12]*pB[0] + pA[13]*pB[1] + pA[14]*pB[2] +
2*pA[15]*pB[3] + 2*pA[16]*pB[4] + 2*pA[17]*pB[5]) +
2*pB[3]*(pA[18]*pB[0] + pA[19]*pB[1] + pA[20]*pB[2] +
2*pA[21]*pB[3] + 2*pA[22]*pB[4] + 2*pA[23]*pB[5]) +
2*pB[4]*(pA[24]*pB[0] + pA[25]*pB[1] + pA[26]*pB[2] +
2*pA[27]*pB[3] + 2*pA[28]*pB[4] + 2*pA[29]*pB[5]) +
2*pB[5]*(pA[30]*pB[0] + pA[31]*pB[1] + pA[32]*pB[2] +
2*pA[33]*pB[3] + 2*pA[34]*pB[4] + 2*pA[35]*pB[5]);
	else if (fNumSD == 1)
		return fArray[0]*fArray[0]*A[0];
	else
		throw eGeneralFail;
		
	/* dummy */
	return 0;
}

/* outer product */
void dSymMatrixT::Outer(const dArrayT& v)
{
#if __option (extended_errorcheck)
	/* dimension checks */
	if (fNumSD != v.Length()) throw eSizeMismatch;
#endif

	double *pthis = fArray;
	double    *pv = v.Pointer();
	if (fNumSD == 2)
	{
		*pthis++ = pv[0]*pv[0];
		*pthis++ = pv[1]*pv[1];
		*pthis   = pv[0]*pv[1];
	}
	else if (fNumSD == 3)
	{	
		*pthis++ = pv[0]*pv[0];
		*pthis++ = pv[1]*pv[1];
		*pthis++ = pv[2]*pv[2];
		
		*pthis++ = pv[1]*pv[2];
		*pthis++ = pv[0]*pv[2];
		*pthis   = pv[0]*pv[1];
	}
	else
		*pthis = pv[0]*pv[0];
}

/* matrix-matrix multiplication */
void dSymMatrixT::MultAB(const dSymMatrixT& A, const dSymMatrixT& B)
{
#if __option (extended_errorcheck)
	/* dimension checks */
	if (  fNumSD != A.fNumSD ||
	    A.fNumSD != B.fNumSD) throw eSizeMismatch;
#endif

	if (fNumSD == 2)
	{
		double* pA = A.fArray;
		double* pB = B.fArray;
	
		fArray[0] = pA[0]*pB[0] + pA[2]*pB[2];
		fArray[1] = pA[1]*pB[1] + pA[2]*pB[2];
		fArray[2] = pA[2]*pB[1] + pA[0]*pB[2];
	}
	else if (fNumSD == 3)
	{
		double* pA = A.fArray;
		double* pB = B.fArray;

		fArray[0] = pA[0]*pB[0] + pA[4]*pB[4] + pA[5]*pB[5];
		fArray[1] = pA[1]*pB[1] + pA[3]*pB[3] + pA[5]*pB[5];
		fArray[2] = pA[2]*pB[2] + pA[3]*pB[3] + pA[4]*pB[4];
		fArray[3] = pA[3]*pB[2] + pA[1]*pB[3] + pA[5]*pB[4];
		fArray[4] = pA[4]*pB[2] + pA[5]*pB[3] + pA[0]*pB[4];
		fArray[5] = pA[5]*pB[1] + pA[4]*pB[3] + pA[0]*pB[5];
	}
	else if (fNumSD == 1)
		fArray[0] = A[0]*B[0];
	else
		throw eGeneralFail;
}

void dSymMatrixT::MultATA(const dMatrixT& A)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (A.Rows() != A.Cols() ||
	    A.Rows() != fNumSD) throw eSizeMismatch;
#endif

	double* pA = A.Pointer();

	if (fNumSD == 2)
	{
		fArray[0] = pA[0]*pA[0] + pA[1]*pA[1];
		fArray[1] = pA[2]*pA[2] + pA[3]*pA[3];
		fArray[2] = pA[0]*pA[2] + pA[1]*pA[3];
	}
	else if (fNumSD == 3)
	{
		fArray[0] = pA[0]*pA[0] + pA[1]*pA[1] + pA[2]*pA[2];
		fArray[1] = pA[3]*pA[3] + pA[4]*pA[4] + pA[5]*pA[5];
		fArray[2] = pA[6]*pA[6] + pA[7]*pA[7] + pA[8]*pA[8];

		fArray[3] = pA[3]*pA[6] + pA[4]*pA[7] + pA[5]*pA[8];
		fArray[4] = pA[0]*pA[6] + pA[1]*pA[7] + pA[2]*pA[8];
		fArray[5] = pA[0]*pA[3] + pA[1]*pA[4] + pA[2]*pA[5];
	}
	else if (fNumSD == 1)
		fArray[0] = pA[0]*pA[0];
	else
		throw eGeneralFail;
}

void dSymMatrixT::MultAAT(const dMatrixT& A)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (A.Rows() != A.Cols() ||
	    A.Rows() != fNumSD) throw eSizeMismatch;
#endif

	double* pA = A.Pointer();

	if (fNumSD == 2)
	{
		fArray[0] = pA[0]*pA[0] + pA[2]*pA[2];
		fArray[1] = pA[1]*pA[1] + pA[3]*pA[3];
		fArray[2] = pA[0]*pA[1] + pA[2]*pA[3];
	}
	else if (fNumSD == 3)
	{
		fArray[0] = pA[0]*pA[0] + pA[3]*pA[3] + pA[6]*pA[6];
		fArray[1] = pA[1]*pA[1] + pA[4]*pA[4] + pA[7]*pA[7];
		fArray[2] = pA[2]*pA[2] + pA[5]*pA[5] + pA[8]*pA[8];

		fArray[3] = pA[1]*pA[2] + pA[4]*pA[5] + pA[7]*pA[8];
		fArray[4] = pA[0]*pA[2] + pA[3]*pA[5] + pA[6]*pA[8];
		fArray[5] = pA[0]*pA[1] + pA[3]*pA[4] + pA[6]*pA[7];
	}
	else if (fNumSD == 1)
		fArray[0] = pA[0]*pA[0];
	else
		throw eGeneralFail;
}

/* matrix-matrix-matrix operations, ie. tensor basis transformations */
void dSymMatrixT::MultQBQT(const dMatrixT& Q, const dSymMatrixT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (  fNumSD != Q.Rows() ||
	    B.fNumSD != Q.Cols() ||
	    B.fNumSD != fNumSD) throw eSizeMismatch;
#endif

	double* b = B.Pointer();
	double* q = Q.Pointer();

	if (fNumSD == 2)
	{
		fArray[0] = b[0]*q[0]*q[0] + 2.0*b[2]*q[0]*q[2] + b[1]*q[2]*q[2];
		fArray[1] = b[0]*q[1]*q[1] + 2.0*b[2]*q[1]*q[3] + b[1]*q[3]*q[3];
		fArray[2] = b[0]*q[0]*q[1] + b[2]*q[1]*q[2] +
		            b[2]*q[0]*q[3] + b[1]*q[2]*q[3];
	}
	else if (fNumSD == 3)
	{	
		fArray[0] = b[0]*q[0]*q[0] + 2*b[5]*q[0]*q[3] + b[1]*q[3]*q[3] +
		          2*b[4]*q[0]*q[6] + 2*b[3]*q[3]*q[6] + b[2]*q[6]*q[6];
		fArray[1] = b[0]*q[1]*q[1] + 2*b[5]*q[1]*q[4] + b[1]*q[4]*q[4] +
		          2*b[4]*q[1]*q[7] + 2*b[3]*q[4]*q[7] + b[2]*q[7]*q[7];	
		fArray[2] = b[0]*q[2]*q[2] + 2*b[5]*q[2]*q[5] + b[1]*q[5]*q[5] +
		          2*b[4]*q[2]*q[8] + 2*b[3]*q[5]*q[8] + b[2]*q[8]*q[8];
		
		fArray[3] = q[5]*(b[5]*q[1] + b[1]*q[4] + b[3]*q[7]) +
		            q[2]*(b[0]*q[1] + b[5]*q[4] + b[4]*q[7]) +
		            q[8]*(b[4]*q[1] + b[3]*q[4] + b[2]*q[7]);
		fArray[4] = q[5]*(b[5]*q[0] + b[1]*q[3] + b[3]*q[6]) +
		            q[2]*(b[0]*q[0] + b[5]*q[3] + b[4]*q[6]) +
		            q[8]*(b[4]*q[0] + b[3]*q[3] + b[2]*q[6]);
		fArray[5] = q[4]*(b[5]*q[0] + b[1]*q[3] + b[3]*q[6]) +
		            q[1]*(b[0]*q[0] + b[5]*q[3] + b[4]*q[6]) +
		            q[7]*(b[4]*q[0] + b[3]*q[3] + b[2]*q[6]);	
}	
	else if (fNumSD == 1)
		fArray[0] = b[0]*q[0]*q[0];
	else
		throw eGeneralFail;
}

void dSymMatrixT::MultQTBQ(const dMatrixT& Q, const dSymMatrixT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (  fNumSD != Q.Cols() ||
	    B.fNumSD != Q.Rows() ||
	    B.fNumSD != fNumSD) throw eSizeMismatch;
#endif

	double* b = B.Pointer();
	double* q = Q.Pointer();

	if (fNumSD == 2)
	{
		fArray[0] = b[0]*q[0]*q[0] + 2.0*b[2]*q[0]*q[1] + b[1]*q[1]*q[1];
		fArray[1] = b[0]*q[2]*q[2] + 2.0*b[2]*q[2]*q[3] + b[1]*q[3]*q[3];
		fArray[2] = b[0]*q[0]*q[2] + b[2]*q[1]*q[2] +
		            b[2]*q[0]*q[3] + b[1]*q[1]*q[3];
	}
	else if (fNumSD == 3)
	{
		fArray[0] = b[0]*q[0]*q[0] + 2*b[5]*q[0]*q[1] + b[1]*q[1]*q[1] +
		          2*b[4]*q[0]*q[2] + 2*b[3]*q[1]*q[2] + b[2]*q[2]*q[2];
		fArray[1] = b[0]*q[3]*q[3] + 2*b[5]*q[3]*q[4] + b[1]*q[4]*q[4] +
		          2*b[4]*q[3]*q[5] + 2*b[3]*q[4]*q[5] + b[2]*q[5]*q[5];
		fArray[2] = b[0]*q[6]*q[6] + 2*b[5]*q[6]*q[7] + b[1]*q[7]*q[7] +
		          2*b[4]*q[6]*q[8] + 2*b[3]*q[7]*q[8] + b[2]*q[8]*q[8];

		fArray[3] = q[5]*(b[4]*q[6] + b[3]*q[7] + b[2]*q[8]) +
		            q[4]*(b[5]*q[6] + b[1]*q[7] + b[3]*q[8]) +
		            q[3]*(b[0]*q[6] + b[5]*q[7] + b[4]*q[8]);
		fArray[4] = (b[0]*q[0] + b[5]*q[1] + b[4]*q[2])*q[6] +
		            (b[5]*q[0] + b[1]*q[1] + b[3]*q[2])*q[7] +
		            (b[4]*q[0] + b[3]*q[1] + b[2]*q[2])*q[8];
		fArray[5] = (b[0]*q[0] + b[5]*q[1] + b[4]*q[2])*q[3] +
		            (b[5]*q[0] + b[1]*q[1] + b[3]*q[2])*q[4] +
		            (b[4]*q[0] + b[3]*q[1] + b[2]*q[2])*q[5];
	}
	else if (fNumSD == 1)
		fArray[0] = b[0]*q[0]*q[0];
	else
		throw eGeneralFail;
}

/* matrix-vector multiplication */
void dSymMatrixT::Multx(const dArrayT& x, dArrayT& b) const
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (x.Length() != b.Length() ||
	    b.Length() != fNumSD) throw eSizeMismatch;
#endif

	double* px = x.Pointer();

	if (fNumSD == 2)
	{
		b[0] = fArray[0]*px[0] + fArray[2]*px[1];
		b[1] = fArray[2]*px[0] + fArray[1]*px[1];		
	}
	else if (fNumSD == 3)
	{
		b[0] = fArray[0]*px[0] + fArray[5]*px[1] + fArray[4]*px[2];
		b[1] = fArray[5]*px[0] + fArray[1]*px[1] + fArray[3]*px[2];		
		b[2] = fArray[4]*px[0] + fArray[3]*px[1] + fArray[2]*px[2];		
	}
	else if (fNumSD == 1)
		b[0] = fArray[0]*px[0];
	else
		throw eGeneralFail;
}

/* vector-matrix-vector product */
double dSymMatrixT::MultmBn(const dArrayT& m, const dArrayT& n) const
{
/* dimension check */
#if __option (extended_errorcheck)
	if (m.Length() != fNumSD ||
	    n.Length() != fNumSD) throw eSizeMismatch;
#endif

	double *pn = n.Pointer();
	double *pm = m.Pointer();

	if (fNumSD == 2)
		return (fArray[0]*pm[0] + fArray[2]*pm[1])*pn[0] +
		       (fArray[2]*pm[0] + fArray[1]*pm[1])*pn[1];
	else if (fNumSD == 3)
		return (fArray[0]*pm[0] + fArray[5]*pm[1] + fArray[4]*pm[2])*pn[0] +
(fArray[5]*pm[0] + fArray[1]*pm[1] + fArray[3]*pm[2])*pn[1] +
(fArray[4]*pm[0] + fArray[3]*pm[1] + fArray[2]*pm[2])*pn[2];
	else if (fNumSD == 1)
		return fArray[0]*pm[0]*pn[0];
	else
		throw eGeneralFail;

	return 0.0;
}

/**************************************************************************
* Private
**************************************************************************/

int dSymMatrixT::Eigenvalues3D(dArrayT& evals, bool sort_descending,
	int max_iterations) const
{
#if __option(extended_errorcheck)
	if (Rows() != evals.Length() ||
evals.Length() != 3) throw eSizeMismatch;
#endif

	double a[3], b[3], z[3];
	double* d = evals.Pointer();

	/* init 1D arrays */
	a[0] = (*this)[5]; // 1,2
	a[1] = (*this)[3]; // 2,3
	a[2] = (*this)[4]; // 1,3

	b[0] = d[0] = (*this)[0]; // 1,1
	b[1] = d[1] = (*this)[1]; // 2,2
	b[2] = d[2] = (*this)[2]; // 3,3

	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = 0.0;

	int permute[4] = {1, 2, 0, 1};
	int rotations = 0;
	int iterations = 0;
	int min_iterations = 4;
	double small = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);
	while (iterations++ < 50 && small != 0.0)
	{
		double threshold = (iterations < min_iterations) ? 0.011*small : 0.0;

		/* rotations */
		for (int i = 0; i < 3; i++)
		{	
			int j = permute[i];
			int k = permute[j];
	
			double aij = a[i];
			double g = 100.0*fabs(aij);
			if (fabs(d[i]) + g != fabs(d[i]) || fabs(d[j]) + g != fabs(d[j]))
			{
				if (fabs(aij) > threshold)
				{
					a[i] = 0.0;
					double h = d[j] - d[i];
					double t;
					if (fabs(h) + g == fabs(h))
						t = aij/h;
					else
					{
						double ha = h/aij;
						t = d_sign(2.0,ha)/(fabs(ha) + sqrt(4.0 + ha*ha));
					}
					
					/* set rotation */
					double c = 1.0/sqrt(1.0 + t*t);
					double s = t*c;
					double tau = s/(1.0 + c);
					
					/* rotate diagonals */
					h = t*aij;
					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;
					
					/* rotate off-diagonals */
					h = a[j];
					g = a[k];
					a[j] = h + s*(g - h*tau);
					a[k] = g - s*(h + g*tau);
					
					rotations++;
				}			
			}
			else
				a[i] = 0.0;
	
			/* update diagonals */
			double* _b = b;
			double* _d = d;
			double* _z = z;		
			*_b += *_z;
			*_d++ = *_b++;
			*_z++ = 0.0;	
			*_b += *_z;
			*_d++ = *_b++;
			*_z++ = 0.0;	
			*_b += *_z;
			*_d++ = *_b++;
			*_z++ = 0.0;	
		}
		small = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);
	}

	/* check convergence */
	if (iterations >= max_iterations)
	{
		cout << "\n dSymMatrixT::Eigenvalues3D: failed to converge after "
		     << max_iterations << " iterations\n"
		     <<   "     Sum of off-diagonal terms = " << small << endl;
		throw eGeneralFail;
	}
	
	if (sort_descending)
	{
		double tmp;
		if (d[0] < d[1])
		{
			tmp = d[0];
			d[0] = d[1];
			d[1] = tmp;
		}				
		
		if (d[0] < d[2])
		{
			tmp = d[2];
			d[2] = d[1];
			d[1] = d[0];
			d[0] = tmp;
		}
		else if (d[1] < d[2])
		{
			tmp = d[1];
			d[1] = d[2];
			d[2] = tmp;
		}
	}			
	
	return rotations;
}

int dSymMatrixT::Eigensystem3D(dArrayT& evals, dMatrixT& evecs, bool sort_descending,
	int max_iterations) const
{
#if __option(extended_errorcheck)
	if (Rows() != evals.Length() ||
evecs.Rows() != evecs.Cols() ||
evecs.Rows() != evals.Length() ||
evecs.Rows() != 3) throw eSizeMismatch;
#endif

//TEMP
if (sort_descending)
{
	cout << "\n dSymMatrixT::Eigensystem3D: sort_descending not implemented" << endl;
	throw eGeneralFail;
}

	double a[3], b[3], z[3];
	double* d = evals.Pointer();

	/* init 1D arrays */
	a[0] = (*this)[5]; // 1,2
	a[1] = (*this)[3]; // 2,3
	a[2] = (*this)[4]; // 1,3

	b[0] = d[0] = (*this)[0]; // 1,1
	b[1] = d[1] = (*this)[1]; // 2,2
	b[2] = d[2] = (*this)[2]; // 3,3

	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = 0.0;
	
	double* v = evecs.Pointer();
	*v++ = 1.0;
	*v++ = 0.0;
	*v++ = 0.0;
	*v++ = 0.0;
	*v++ = 1.0;
	*v++ = 0.0;
	*v++ = 0.0;
	*v++ = 0.0;
	*v   = 1.0;

	int permute[4] = {1, 2, 0, 1};
	int rotations = 0;
	int iterations = 0;
	int min_iterations = 4;
	double small = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);
	while (iterations++ < 50 && small != 0.0)
	{
		double threshold = (iterations < min_iterations) ? 0.011*small : 0.0;

		/* rotations */
		for (int i = 0; i < 3; i++)
		{	
			int j = permute[i];
			int k = permute[j];
	
			double aij = a[i];
			double g = 100.0*fabs(aij);
			if (fabs(d[i]) + g != fabs(d[i]) || fabs(d[j]) + g != fabs(d[j]))
			{
				if (fabs(aij) > threshold)
				{
					a[i] = 0.0;
					double h = d[j] - d[i];
					double t;
					if (fabs(h) + g == fabs(h))
						t = aij/h;
					else
					{
						double ha = h/aij;
						t = d_sign(2.0,ha)/(fabs(ha) + sqrt(4.0 + ha*ha));
					}
					
					/* set rotation */
					double c = 1.0/sqrt(1.0 + t*t);
					double s = t*c;
					double tau = s/(1.0 + c);
					
					/* rotate diagonals */
					h = t*aij;
					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;
					
					/* rotate off-diagonals */
					h = a[j];
					g = a[k];
					a[j] = h + s*(g - h*tau);
					a[k] = g - s*(h + g*tau);
					
					/* rotate eigenvectors */
					double* v = evecs.Pointer();
					g = v[i];
					h = v[j];
					v[i] = g - s*(h + g*tau);
					v[j] = h + s*(g - h*tau);
					v += 3;
					g = v[i];
					h = v[j];
					v[i] = g - s*(h + g*tau);
					v[j] = h + s*(g - h*tau);
					v += 3;
					g = v[i];
					h = v[j];
					v[i] = g - s*(h + g*tau);
					v[j] = h + s*(g - h*tau);
					
					rotations++;
				}			
			}
			else
				a[i] = 0.0;
	
			/* update diagonals */
			double* _b = b;
			double* _d = d;
			double* _z = z;		
			*_b += *_z;
			*_d++ = *_b++;
			*_z++ = 0.0;	
			*_b += *_z;
			*_d++ = *_b++;
			*_z++ = 0.0;	
			*_b += *_z;
			*_d++ = *_b++;
			*_z++ = 0.0;	
		}

		small = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);
	}

	/* check convergence */
	if (iterations >= max_iterations)
	{
		cout << "\n dSymMatrixT::Eigensystem3D: failed to converge after "
		     << max_iterations << " iterations\n"
		     <<   "     Sum of off-diagonal terms = " << small << endl;
		throw eGeneralFail;
	}

#if 0
	if (sort_descending)
	{
		int map[3];
		double tmp;
		if (d[0] < d[1])
		{
			tmp = d[0];
			d[0] = d[1];
			d[1] = tmp;
		}				
		
		if (d[0] < d[2])
		{
			tmp = d[2];
			d[2] = d[1];
			d[1] = d[0];
			d[0] = tmp;
		}
		else if (d[1] < d[2])
		{
			tmp = d[1];
			d[1] = d[2];
			d[2] = tmp;
		}
	}			
#endif

	return rotations;
}

/* closed form eigenvalues - unstable for repeated roots */
void dSymMatrixT::Eigenvalues3D_Cardano(dArrayT& eigs) const
{
	/* flag */
	bool failed = false;
		
	/* non-const temporary */
	double* tmp = (double*) fArray;
		
	/* apply shift */		
	double shift = (tmp[0] + tmp[1] + tmp[2])/3.0;
	tmp[0] -= shift;
	tmp[1] -= shift;
	tmp[2] -= shift;		
	
	/* compute tensor invariants */
	//double I1 = Trace(); -> becomes zero with the shift
	double I2 =-Invariant2();
	double I3 = Det();
				
	//double b = I2 - I1*I1/3.0;
	//double c =(-2.0*I1*I1*I1 + I1*I2*9.0 - 27.0*I3)/27.0;
	double b = I2;
	double c =-I3;
		
	if (fabs(b) < kSmall)
	{
		/* check c */
		if (c < 0.0)
			failed = true;
		else
		{			
			//double eig = I1/3.0 - pow(c,1.0/3.0);			
			double eig =-pow(c,1.0/3.0);			
			eigs[0] = eig;
			eigs[1] = eig;
			eigs[2] = eig;
		}
	}
	else if (b > 0.0)
		failed = true;
	else
	{
		double m = 2.0*sqrt(-b/3.0);
		double n = 3.0*c/(m*b);	
		if (n > 1.0 && n - 1.0 < 10*kSmall) n = 1.0;
		if (fabs(n) > 1.0)
			failed = true;
		else
		{
			double t = atan2(sqrt(1 - n*n),n)/3.0;
			for (int A = 0; A < 3; A++)
				//eigs[A] = m*cos(t + 2.0*A*Pi/3.0) + I1/3.0;
				eigs[A] = m*cos(t + 2.0*A*Pi/3.0);
		}
	}

	/* failed */
	if (failed)
	{
		eigs[0] = 0.0;
		eigs[1] = 0.0;
		eigs[2] = 0.0;
	}

	/* shift eigenvalues */
	eigs[0] += shift;
	eigs[1] += shift;
	eigs[2] += shift;

	/* restore matrix */
	tmp[0] += shift;
	tmp[1] += shift;
	tmp[2] += shift;		
}
