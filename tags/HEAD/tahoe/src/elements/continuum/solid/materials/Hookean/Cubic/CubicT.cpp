/* $Id: CubicT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "CubicT.h"

#include <iostream.h>

#include "fstreamT.h"
#include "dMatrixT.h"

/* constructor */
CubicT::CubicT(ifstreamT& in, dMatrixT& moduli)
{	
	in >> fC11;			
	in >> fC12;			//add check on the ranges!!!!
	in >> fC44;	
	
	/* set moduli */
	ComputeModuli(moduli, fC11, fC12, fC44);
}

/* I/O operators */
void CubicT::Print(ostream& out) const
{
	out << " C11 . . . . . . . . . . . . . . . . . . . . . . = " << fC11 << '\n';
	out << " C12 . . . . . . . . . . . . . . . . . . . . . . = " << fC12 << '\n';
	out << " C44 . . . . . . . . . . . . . . . . . . . . . . = " << fC44 << '\n';	
}

/*************************************************************************
* Protected
*************************************************************************/

/* print name */
void CubicT::PrintName(ostream& out) const
{
	out << "    Cubic\n";
}

/* compute the symetric Cij reduced index matrix */
void CubicT::ComputeModuli(dMatrixT& moduli,
	double C11, double C12, double C44)
{
	moduli = 0.0;

	if (moduli.Rows() == 3)
	{
		moduli(1,1) = moduli(0,0) = C11;
		moduli(0,1) = C12;
		moduli(2,2) = C44;
	}
	else
	{
		moduli(2,2) = moduli(1,1) = moduli(0,0) = C11;
		moduli(1,2) = moduli(0,1) = moduli(0,2) = C12;
		moduli(5,5) = moduli(4,4) = moduli(3,3) = C44;
	}

	/* symmetric */
	moduli.CopySymmetric();
}
