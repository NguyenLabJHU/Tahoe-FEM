/* $Id: ConstQuadT.cpp,v 1.1 2002-05-31 07:17:12 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "ConstQuadT.h"
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"

/* constructors */
ConstQuadT::ConstQuadT(double A, double B, double C): 
	fA(A),
	fB(B),
	fC(C)
 { }

/* I/O */
void ConstQuadT::Print(ostream& out) const
{
	/* parameters */
        out <<"\n        A = " << fA; 
	out <<"\n        B = " << fB;  
        out <<"\n        C = " << fC  << '\n';
} 

void ConstQuadT::PrintName(ostream& out) const
{
	out << "Constant-Quadratic \n";
}

/* returning values in groups */
dArrayT& ConstQuadT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x;

	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pU++ = Function(x);
	}

	return out;
}

dArrayT& ConstQuadT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x;

	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pU++ = DFunction(x);
	}

	return out;
}

dArrayT& ConstQuadT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x;

	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pU++ = DDFunction(x);
	}
	return out;
}