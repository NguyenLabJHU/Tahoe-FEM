/* $Id: Rotate2DT.cpp,v 1.5 2005-01-13 00:10:32 paklein Exp $ */
/* created: paklein (07/21/1996)                                          */
/* This class provides the functionality to do 2D coordinate              */
/* transformations.                                                       */

#include "Rotate2DT.h"
#include <math.h>
#include "toolboxConstants.h"

/* size parameters */

using namespace Tahoe;

const int kMatrixDim = 2;
const int kRedMatDim = 3;

const double Pi = acos(-1.0);

/*
* constructor
*/
Rotate2DT::Rotate2DT(void): fAngleDeg(0.0), fAngle(0.0), fCost(0.0),
	fSint(0.0), fCos2t(0.0), fSin2t(0.0), fCos4t(0.0), fSin4t(0.0),
	fQ(kMatrixDim), fVec(kMatrixDim), fRedMat(kRedMatDim), fMat(kMatrixDim)
{

}

Rotate2DT::Rotate2DT(double angle): fAngleDeg(angle), fQ(kMatrixDim),
	fVec(kMatrixDim), fRedMat(kRedMatDim), fMat(kMatrixDim)
{
	SetAngle(fAngleDeg);
}

/*
* set the angle field and associated work variables.
*/
void Rotate2DT::SetAngle(double angle)
{
	fAngleDeg = angle;
	fAngle = fAngleDeg*Pi/180.0;

	double temp = fAngle;
	fCost = cos(temp);	
	fSint = sin(temp);
	
	temp *= 2;
	fCos2t = cos(temp);	
	fSin2t = sin(temp);
	
	temp *= 2;
	fCos4t = cos(temp);	
	fSin4t = sin(temp);

	/* rotation matrix */
	fQ(0,0) = fQ(1,1) = fCost;
	fQ(0,1) =-fSint;
	fQ(1,0) = fSint;		
}

/*
* reduced index symmetric matrices
*/
const dSymMatrixT& Rotate2DT::RotateRedMatIn(const dSymMatrixT& redmat)
{
/* dimension checking */
#if __option (extended_errorcheck)
	if (redmat.Length() != kRedMatDim) throw ExceptionT::kGeneralFail;
#endif

	/* Note: This reduced index coordinate transformation requires
	         21 operations, where a full, non-symmetric transformation
	         takes 36 operations, not including the sine's and cosine's,
	         which are computed ahead of time. */

	double sigmean = 0.5*(redmat[0] + redmat[1]);
	double sigdiff = 0.5*(redmat[0] - redmat[1]);
	
	fRedMat[0] = sigmean + sigdiff*fCos2t + redmat[2]*fSin2t;
	fRedMat[1] = sigmean - sigdiff*fCos2t - redmat[2]*fSin2t;
	fRedMat[2] = redmat[2]*fCos2t - sigdiff*fSin2t;
	
	return(fRedMat);
}

const dSymMatrixT& Rotate2DT::RotateRedMatOut(const dSymMatrixT& redmat)
{
/* dimension checking */
#if __option (extended_errorcheck)
	if (redmat.Length() != kRedMatDim) throw ExceptionT::kGeneralFail;
#endif

	/* Note: This reduced index coordinate transformation requires
	         21 operations, where a full, non-symmetric transformation
	         takes 36 operations, not include the sine's and cosine's,
	         which are computed ahead of time. */

	double sigmean = 0.5*(redmat[0] + redmat[1]);
	double sigdiff = 0.5*(redmat[0] - redmat[1]);
	
	fRedMat[0] = sigmean + sigdiff*fCos2t - redmat[2]*fSin2t;
	fRedMat[1] = sigmean - sigdiff*fCos2t + redmat[2]*fSin2t;
	fRedMat[2] = redmat[2]*fCos2t + sigdiff*fSin2t;

	return(fRedMat);
}

/* push an index in */
void Rotate2DT::PushIndexIn(dMatrixT& matrix, int index)
{
/* dimension checking */
#if __option (extended_errorcheck)
	if (matrix.Rows() != kMatrixDim ||
	    matrix.Cols() != kMatrixDim ||
	    index < 0                   ||
	    index >= kMatrixDim) throw ExceptionT::kGeneralFail;
#endif

	/* copy */
	fMat = matrix;

	/* second index */
	if (index == 1)
		matrix.MultAB(fMat,fQ);
	/* first index */
	else
		matrix.MultATB(fQ,fMat);
}

/*
* reduced index 4th order tensors
*/
void Rotate2DT::RotateRedTensorIn(dMatrixT& matrix)
{
/* dimension checking */
#if __option (extended_errorcheck)
	if (matrix.Rows() != kRedMatDim ||
	    matrix.Cols() != kRedMatDim) throw ExceptionT::kGeneralFail;
#endif

	TransformO42D(matrix, 1);
}

void Rotate2DT::RotateRedTensorOut(dMatrixT& matrix)
{
/* dimension checking */
#if __option (extended_errorcheck)
	if (matrix.Rows() != kRedMatDim ||
	    matrix.Cols() != kRedMatDim) throw ExceptionT::kGeneralFail;
#endif

	TransformO42D(matrix,-1);
}

/***********************************************************************
* Private
***********************************************************************/

/*
* Crunch the numbers - transform direction set the direction either
* (+) or (-) fAngle for (1) and (-1), respectively.
*/
void Rotate2DT::TransformO42D(dMatrixT& matrix, int RotateDirection)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;
	double z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, z27, z28;
	double z30, z31, z32, z33, z34, z35, z36, z37, z38, z40, z41;
	double z42, z44, z45, z46;
	//not used
	//double z29, z39, z43;

	/* check */
	if (RotateDirection != 1 &&
	    RotateDirection != -1) throw ExceptionT::kGeneralFail;

	z3 = matrix(0,0);
	z4 = matrix(0,1);
	z5 = matrix(0,2);
	z6 = matrix(1,1);
	z7 = matrix(1,2);
	z8 = matrix(2,2);
	z9 = fCost;
	z10 = fSint*RotateDirection;
	z11 = z10*z10;
	z12 = z11*z10;
	z13 = z11*z11;
	z14 = -2*z4;
	z15 = 6*z4;
	z16 = -4*z8;
	z8 = 4*z8;
	z17 = z9*z9;
	z18 = z9*z17;
	z19 = z17*z17;
	z20 = fCos2t;
	z21 = fCos4t;
	z1 = fSin2t*RotateDirection;
	z2 = fSin4t*RotateDirection;
	z11 = z11*z17;
	z10 = z10*z18;
	z17 = 4*z20;
	z18 = -4*z21;
	z20 = -z21;
	z22 = 4*z21;
	z23 = z21*z8;
	z24 = -z2;
	z25 = z14*z2;
	z26 = z16*z2;
	z27 = z2*z8;
	z28 = z11*z8;
	//z29 = -z3;
	z30 = z13*z3;
	z31 = z19*z3;
	z32 = 2*z1*z3;
	z33 = z2*z3;
	z34 = z20*z3;
	z35 = z24*z3;
	z4 = 2*z4;
	z36 = z21*z4;
	z37 = z2*z4;
	z4 = z11*z4;
	z11 = -4*z5;
	z11 = z10*z11;
	z38 = 4*z5;
	//z39 = z21*z5;
	z40 = z17*z5;
	z41 = z18*z5;
	z5 = z22*z5;
	z42 = z2*z38;
	//z43 = -z6;
	z13 = z13*z6;
	z19 = z19*z6;
	z1 = -2*z1*z6;
	z44 = z2*z6;
	z20 = z20*z6;
	z24 = z24*z6;
	z45 = -4*z7;
	z2 = z2*z45;
	z46 = 4*z7;
	z10 = z10*z46;
	z21 = z21*z7;
	z17 = z17*z7;
	z18 = z18*z7;
	z7 = z22*z7;
	z9 = z12*z9;
	z12 = z38*z9;
	z9 = z45*z9;
	z4 = z28 + z4;
	z2 = z2 + z20 + z23 + z3 + z34 + z36 + z42 + z6;
	z1 = z1 + z17 + z32 + z40;
	z3 = z10 + z12 + z19 + z30 + z4;
	z4 = z11 + z13 + z31 + z4 + z9;
	z6 = z15 + z16 + z2;
	z2 = z14 + z2 + z8;
	z5 = z1 + z18 + z25 + z26 + z33 + z44 + z5;
	z1 = z1 + z24 + z27 + z35 + z37 + z41 + z7;

	matrix(0,0) = z4;
	matrix(0,2) = matrix(2,0) = z5/8;
	matrix(2,2) = z2/8;
	matrix(0,1) = matrix(1,0) = z6/8;
	matrix(1,2) = matrix(2,1) = z1/8;
	matrix(1,1) = z3;	
}
