/* $Id: CubicT.cpp,v 1.6 2004-06-28 22:41:22 hspark Exp $ */
/* created: paklein (06/11/1997) */

#include "CubicT.h"

#include <iostream.h>

#include "ifstreamT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* constructor */
CubicT::CubicT(ifstreamT& in)
{	
	in >> fC11;			
	in >> fC12;			//add check on the ranges!!!!
	in >> fC44;	
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
void CubicT::ComputeModuli(dMatrixT& moduli)
{
	if (moduli.Rows() == 6)
	{
		moduli = 0.0;
		moduli(2,2) = moduli(1,1) = moduli(0,0) = fC11;
		moduli(1,2) = moduli(0,1) = moduli(0,2) = fC12;
		moduli(5,5) = moduli(4,4) = moduli(3,3) = fC44;

		/* symmetric */
		moduli.CopySymmetric();
	}
	else
	{
		cout << "\n CubicT::ComputeModuli: 3D only" << endl;
		throw ExceptionT::kSizeMismatch;
	}
}

void CubicT::ComputeModuli2D(dMatrixT& moduli, Material2DT::ConstraintOptionT constraint) const
{
	if (moduli.Rows() == 3)
	{
		/* reset moduli for plane stress */
		double C11, C12;
		if (constraint == Material2DT::kPlaneStress)
		{
			C11 = fC11 - (fC12*fC12/fC11);
			C12 = fC12*(fC11 - fC12)/fC11;
		}
		else
		{
			C11 = fC11;
			C12 = fC12;
		}
		
		moduli = 0.0;
		moduli(1,1) = moduli(0,0) = C11;
		moduli(0,1) = moduli(1,0) = C12;
		moduli(2,2) = fC44;
	}
	else throw ExceptionT::kSizeMismatch;
}

/* scale factor for constrained dilatation */
double CubicT::DilatationFactor2D(Material2DT::ConstraintOptionT constraint) const
{
	/* scale thermal strain */
	if (constraint == Material2DT::kPlaneStrain)
		return 1.0 + (fC12/(fC12 + fC11));
	else
		 return 1.0;
}
