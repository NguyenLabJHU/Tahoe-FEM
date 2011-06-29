/* $Id: BoyceBaseT.cpp,v 1.1 2007-03-08 18:13:01 tdnguye Exp $ */
/* created: TDN (01/22/2000) */
#include "BoyceBaseT.h"

using namespace Tahoe;

/* constructor */
BoyceBaseT::BoyceBaseT(void):
	ParameterInterfaceT("Boyce_BaseT"),
	fSpectralDecompRef(3)
{
}

/*initializes history variable */
void  BoyceBaseT::PointInitialize(void)
{
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
		      
			  fFv.Identity();
			  fFv_n.Identity();

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void BoyceBaseT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fFv_n = fFv;

		/* write to storage */
		Store(element, ip);
	}
}

void BoyceBaseT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fFv = fFv_n;
		
		/* write to storage */
		Store(element, ip);
	}
}

/* form of tangent matrix */
GlobalT::SystemTypeT BoyceBaseT::TangentType(void) const
{
	/* symmetric by default */
	return GlobalT::kNonSymmetric;
}

void BoyceBaseT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}
void BoyceBaseT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;

}

/*************************************************************************
 * Private
 *************************************************************************/

/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void BoyceBaseT::MixedRank4_2D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 2 ||
	    b.Length() != 2 ||
	    rank4_ab.Rows() != 3 ||
	    rank4_ab.Cols() != 3) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = B(1.);
//	z4 = B(2.);

	z1 = a[0];
	z2 = a[1];
	z3 = b[0];
	z4 = b[1];

	z5 = z1*z1;
	z6 = z2*z2;
	z7 = z3*z3;
	z8 = 2.*z1*z2*z3*z4;
	z9 = z4*z4;
	z3 = 2.*z3*z4;
	z4 = z3*z5;
	z3 = z3*z6;
	z10 = 2.*z1*z2*z7;
	z11 = 2.*z5*z7;
	z7 = z6*z7;
	z1 = 2.*z1*z2*z9;
	z2 = z5*z9;
	z5 = 2.*z6*z9;
	z4 = z10 + z4;
	z1 = z1 + z3;
	z2 = z2 + z7 + z8;
	z3 = 0.5*z4;
	z1 = 0.5*z1;
	z2 = 0.5*z2;

	//{{z11, z8, z3}, 
	// {z8, z5, z1}, 
	// {z3, z1, z2}}

	double* p = rank4_ab.Pointer();
	*p++ = z11;
    *p++ = z8;
    *p++ = z3;
    *p++ = z8;
    *p++ = z5;
    *p++ = z1;
    *p++ = z3;
    *p++ = z1;
    *p   = z2;
}

void BoyceBaseT::MixedRank4_3D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 3 ||
	    b.Length() != 3 ||
	    rank4_ab.Rows() != 6 ||
	    rank4_ab.Cols() != 6) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = A(3.);
//	z4 = B(1.);
//	z5 = B(2.);
//	z6 = B(3.);

	z1 = a[0];	
	z2 = a[1];
	z3 = a[2];
	z4 = b[0];
	z5 = b[1];
	z6 = b[2];
	z7 = z1*z1;
	z8 = z2*z2;
	z9 = z3*z3;
	z10 = 2.*z1*z4;
	z11 = z10*z2;
	z12 = z2*z4;
	z13 = z1*z12;
	z14 = z4*z4;
	z15 = z10*z5;
	z15 = z15*z3;
	z16 = z11*z5;
	z17 = z12*z5;
	z18 = 2.*z17;
	z17 = z17*z3;
	z18 = z18*z3;
	z19 = z1*z4*z5;
	z19 = z19*z3;
	z20 = z5*z5;
	z11 = z11*z6;
	z13 = z13*z6;
	z10 = z10*z3*z6;
	z12 = z12*z3*z6;
	z21 = 2.*z12;
	z22 = z1*z2*z5*z6;
	z23 = 2.*z22;
	z24 = z1*z3*z5*z6;
	z25 = 2.*z24;
	z26 = 2.*z2*z3*z5*z6;
	z27 = z6*z6;
	z28 = 2.*z1*z14;
	z29 = 2.*z1*z2;
	z30 = z14*z2;
	z11 = z11 + z15;
	z15 = z1*z20*z3;
	z31 = 2.*z2*z20*z3;
	z18 = z18 + z23;
	z21 = z21 + z25;
	z23 = z1*z2*z27;
	z1 = 2.*z1*z27*z3;
	z25 = 2.*z2*z27*z3;
	z2 = z2*z28;
	z28 = z28*z3;
	z29 = z20*z29;
	z3 = z3*z30;
	z30 = 2.*z14*z7;
	z32 = z20*z7;
	z33 = z27*z7;
	z34 = 2.*z4*z5*z7;
	z35 = 2.*z4*z6*z7;
	z7 = z5*z6*z7;
	z36 = z14*z8;
	z37 = 2.*z20*z8;
	z38 = z27*z8;
	z39 = 2.*z4*z5*z8;
	z40 = z4*z6*z8;
	z8 = 2.*z5*z6*z8;
	z14 = z14*z9;
	z20 = z20*z9;
	z27 = 2.*z27*z9;
	z41 = z4*z5*z9;
	z4 = 2.*z4*z6*z9;
	z5 = 2.*z5*z6*z9;
	z6 = 0.5*z11;
	z9 = 0.5*z18;
	z11 = 0.5*z21;
	z2 = z2 + z34;
	z18 = z28 + z35;
	z3 = z13 + z19 + z3 + z7;
	z7 = z16 + z32 + z36;
	z13 = z29 + z39;
	z15 = z15 + z17 + z22 + z40;
	z8 = z31 + z8;
	z14 = z10 + z14 + z33;
	z17 = z20 + z26 + z38;
	z12 = z12 + z23 + z24 + z41;
	z1 = z1 + z4;
	z4 = z25 + z5;
	z2 = 0.5*z2;
	z5 = 0.5*z18;
	z3 = 0.5*z3;
	z7 = 0.5*z7;
	z13 = 0.5*z13;
	z15 = 0.5*z15;
	z8 = 0.5*z8;
	z14 = 0.5*z14;
	z17 = 0.5*z17;
	z12 = 0.5*z12;
	z1 = 0.5*z1;
	z4 = 0.5*z4;
	
	//{{z30, z16, z10,  z6,  z5,  z2}, 
	// {z16, z37, z26,  z8,  z9, z13}, 
	// {z10, z26, z27,  z4,  z1, z11}, 
	// { z6,  z8,  z4, z17, z12, z15}, 
	// { z5,  z9,  z1, z12, z14,  z3},
	// { z2, z13, z11, z15,  z3,  z7}}
	
	double* p = rank4_ab.Pointer();
    *p++ = z30;
    *p++ = z16;
    *p++ = z10;
    *p++ = z6;
    *p++ = z5;
    *p++ = z2;
    *p++ = z16;
    *p++ = z37;
    *p++ = z26;
    *p++ = z8;
    *p++ = z9;
    *p++ = z13;
    *p++ = z10;
    *p++ = z26;
    *p++ = z27;
    *p++ = z4;
    *p++ = z1;
    *p++ = z11;
    *p++ = z6;
    *p++ = z8;
    *p++ = z4;
    *p++ = z17;
    *p++ = z12;
    *p++ = z15;
    *p++ = z5;
    *p++ = z9;
    *p++ = z1;
    *p++ = z12;
    *p++ = z14;
    *p++ = z3;
    *p++ = z2;
    *p++ = z13;
    *p++ = z11;
    *p++ = z15;
    *p++ = z3;
    *p  = z7;
}

/* accept parameter list */
void BoyceBaseT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);
}

