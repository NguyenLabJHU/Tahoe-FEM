/* $Id: TensorTransformT.cpp,v 1.6 2004-06-26 06:04:50 paklein Exp $ */
/* created: paklein (07/02/1996) */
#include "TensorTransformT.h"

using namespace Tahoe;

/* constructor */
TensorTransformT::TensorTransformT(int dim):
	fRank2(dim),
	fRank4(dSymMatrixT::NumValues(dim)),
	fPull(dim),
	fTransform(dSymMatrixT::NumValues(dim)),
	fOuter(dim),
	fRedMat(dSymMatrixT::NumValues(dim))
{	

}

/* rank 4 tensor transformations */
const dMatrixT& TensorTransformT::PushForward(const dMatrixT& fwd, const dMatrixT& a)
{
	/* copy to internal */
	fRank4 = a;
	
	/* coordinate transformation */
	if (fwd.Rows() == 2)
		FFFFC_2D_Z(fwd, fRank4);
	else if (fwd.Rows() == 3)
		FFFFC_3D(fwd, fRank4);
	else if (fwd.Rows() == 1) {
		double F = fwd[0];
		fRank4[0] *= F*F*F*F;
	}
	else
		ExceptionT::GeneralFail("TensorTransformT::PushForward");
		
	return fRank4;
}

const dMatrixT& TensorTransformT::PullBack(const dMatrixT& fwd, const dMatrixT& a)
{
	/* inverse mapping */
	fPull.Inverse(fwd);

	/* copy to internal */
	fRank4 = a;

	/* coordinate transformation */
	if (fwd.Rows() == 2)
		FFFFC_2D_Z(fPull, fRank4);
	else if (fwd.Rows() == 3)
		FFFFC_3D(fPull, fRank4);
	else if (fwd.Rows() == 1) {
		double F = fPull[0];
		fRank4[0] *= F*F*F*F;
	}
	else
		ExceptionT::GeneralFail("TensorTransformT::PullBack");
			
	return fRank4;
}

/***********************************************************************
* Protected
***********************************************************************/

void TensorTransformT::FFFFC_2D(const dMatrixT& F, dMatrixT& C)
{
#if __option(extended_errorcheck)
	/* consistency */
	if (C.Rows() != dSymMatrixT::NumValues(F.Rows()))
		throw ExceptionT::kSizeMismatch;
#endif

	int nsd = F.Rows();
	dSymMatrixT	coltemp;

	/* compute tranformation matrix */
	dArrayT Fi1(nsd,F(0));
	dArrayT Fi2(nsd,F(1));

	coltemp.Set(nsd,fTransform(0));
	fOuter.Outer(Fi1,Fi1);
	coltemp.FromMatrix(fOuter);

	coltemp.Set(nsd,fTransform(1));
	fOuter.Outer(Fi2,Fi2);
	coltemp.FromMatrix(fOuter);

	coltemp.Set(nsd,fTransform(2));
	fOuter.Outer(Fi2,Fi1);
	fOuter.Symmetrize();
	coltemp.FromMatrix(fOuter);
	coltemp *= 2.0;
	
	/* compute transformed tensor */
	fRedMat.MultQBQT(fTransform, C);
}

void TensorTransformT::FFFFC_3D(const dMatrixT& F, dMatrixT& C)
{
#if __option(extended_errorcheck)
	/* consistency */
	if (C.Rows() != dSymMatrixT::NumValues(F.Rows()))
		throw ExceptionT::kSizeMismatch;
#endif

	int nsd = F.Rows();
	dSymMatrixT	coltemp;

	/* compute tranformation matrix */
	dArrayT Fi1(nsd,F(0));
	dArrayT Fi2(nsd,F(1));
	dArrayT Fi3(nsd,F(2));

	coltemp.Set(nsd,fTransform(0));
	fOuter.Outer(Fi1,Fi1);
	coltemp.FromMatrix(fOuter);

	coltemp.Set(nsd,fTransform(1));
	fOuter.Outer(Fi2,Fi2);
	coltemp.FromMatrix(fOuter);

	coltemp.Set(nsd,fTransform(2));
	fOuter.Outer(Fi3,Fi3);
	coltemp.FromMatrix(fOuter);

	coltemp.Set(nsd,fTransform(3));
	fOuter.Outer(Fi2,Fi3);
	fOuter.Symmetrize();
	coltemp.FromMatrix(fOuter);
	coltemp *= 2.0;

	coltemp.Set(nsd,fTransform(4));
	fOuter.Outer(Fi1,Fi3);
	fOuter.Symmetrize();
	coltemp.FromMatrix(fOuter);
	coltemp *= 2.0;

	coltemp.Set(nsd,fTransform(5));
	fOuter.Outer(Fi1,Fi2);
	fOuter.Symmetrize();
	coltemp.FromMatrix(fOuter);
	coltemp *= 2.0;
	
	/* compute transformed tensor */
	fRedMat.MultQBQT(fTransform, C);
	C = fRedMat;
}

void TensorTransformT::FFFFC_2D_Z(const dMatrixT& F, dMatrixT& C) const
{
#if __option(extended_errorcheck)
	/* consistency */
	if (C.Rows() != dSymMatrixT::NumValues(F.Rows()))
		throw ExceptionT::kSizeMismatch;
#endif

	double	z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;
	double	z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, z27, z28, z29, z30;
	double	z31, z32, z33, z34, z35, z36, z37, z38, z39, z40, z41, z42, z43, z44, z45;
	double	z46, z47, z48, z49, z50, z51;

	z1 = C(0,0);
	z2 = C(0,1);
	z3 = C(0,2);
	z4 = C(1,1);
	z5 = C(1,2);
	z6 = C(2,2);
	z7 = F(0,0);

	z8 = F(0,1);
	z9 = F(1,0);
	z10 = F(1,1);
	z11 = pow(z10,2);
	z12 = pow(z10,3);
	z13 = pow(z10,4);
	z14 = pow(z7,2);
	z15 = pow(z7,3);
	z16 = pow(z7,4);
	z17 = pow(z8,2);
	z18 = pow(z8,3);
	z19 = pow(z8,4);
	z20 = z10*z7*z8*z9;
	z21 = 2*z20;
	z22 = z2*z21;
	z20 = z20*z6;
	z20 = 4*z20;
	z21 = z21*z6;
	z23 = pow(z9,2);
	z24 = pow(z9,3);
	z25 = pow(z9,4);
	z16 = z1*z16;
	z26 = z14*z2;
	z27 = z1*z14*z23;
	z28 = z2*z23;
	z25 = z1*z25;
	z29 = z10*z3;
	z30 = z11*z26;
	z31 = 2*z17*z26;
	z32 = 2*z11*z28;
	z33 = z17*z28;
	z34 = z15*z29;
	z35 = 4*z24*z29;
	z13 = z13*z4;
	z36 = z11*z17*z4;
	z37 = z10*z18*z4;
	z19 = z19*z4;
	z38 = z14*z6;
	z39 = z11*z38;
	z40 = 4*z17*z38;
	z41 = z23*z6;
	z42 = 4*z11*z41;
	z43 = z17*z41;
	z44 = z1*z24*z7;
	z45 = 3*z23*z29*z7;
	z46 = z5*z7;
	z47 = z12*z46;
	z48 = 3*z10*z17*z46;
	z49 = 4*z18*z46;
	z3 = z3*z8;
	z50 = 4*z15*z3;
	z24 = z24*z3;
	z23 = 2*z23*z3*z7;
	z51 = z10*z8;
	z26 = z26*z51;
	z28 = z28*z51;
	z51 = 2*z51;
	z38 = z38*z51;
	z41 = z41*z51;
	z4 = z12*z4*z8;
	z46 = 2*z11*z46*z8;
	z1 = z1*z15*z9;
	z15 = 2*z14*z29*z9;
	z5 = z5*z9;
	z12 = 4*z12*z5;
	z10 = 2*z10*z17*z5;
	z18 = z18*z5;
	z5 = 3*z11*z5*z8;
	z3 = 3*z14*z3*z9;
	z7 = z7*z9;
	z8 = z11*z7;
	z9 = z17*z7;
	z7 = z2*z7;
	z11 = z2*z8;
	z8 = 2*z6*z8;
	z2 = z2*z9;
	z6 = 2*z6*z9;
	z9 = z16 + z19 + z31 + z40 + z49 + z50;
	z12 = z12 + z13 + z25 + z32 + z35 + z42;
	z10 = z10 + z15 + z23 + z27 + z36 + z46;
	z4 = z11 + z24 + z28 + z4 + z41 + z44 + z45 + z47 + z5 + z8;
	z1 = z1 + z18 + z2 + z26 + z3 + z34 + z37 + z38 + z48 + z6;
	z2 = z10 + z20 + z30 + z33;
	z3 = z10 + z21 + z22 + z39 + z43;

	C(0,0) = z9;
	C(2,2) = z3;
	C(1,1) = z12;
	C(0,2) = C(2,0) = z1;
	C(0,1) = C(1,0) = z2;
	C(1,2) = C(2,1) = z4;
}
