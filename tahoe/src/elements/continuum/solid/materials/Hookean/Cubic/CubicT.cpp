/* $Id: CubicT.cpp,v 1.4.34.2 2004-03-02 17:46:13 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "CubicT.h"

#include <iostream.h>

#include "fstreamT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* constructor */
CubicT::CubicT(ifstreamT& in):
	ParameterInterfaceT("cubic")
{	
	in >> fC11;			
	in >> fC12;			//add check on the ranges!!!!
	in >> fC44;	
}

CubicT::CubicT(void): 
	ParameterInterfaceT("cubic"),
	fC11(0.0), fC12(0.0), fC44(0.0)
{

}

/* I/O operators */
void CubicT::Print(ostream& out) const
{
	out << " C11 . . . . . . . . . . . . . . . . . . . . . . = " << fC11 << '\n';
	out << " C12 . . . . . . . . . . . . . . . . . . . . . . = " << fC12 << '\n';
	out << " C44 . . . . . . . . . . . . . . . . . . . . . . = " << fC44 << '\n';	
}

/* describe the parameters needed by the interface */
void CubicT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	list.AddParameter(fC11, "C11");
	list.AddParameter(fC12, "C12");
	list.AddParameter(fC44, "C44");
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

void CubicT::ComputeModuli2D(dMatrixT& moduli, SolidMaterialT::ConstraintT constraint) const
{
	if (moduli.Rows() == 3)
	{
		/* reset moduli for plane stress */
		double C11, C12;
		if (constraint == SolidMaterialT::kPlaneStress)
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
double CubicT::DilatationFactor2D(SolidMaterialT::ConstraintT constraint) const
{
	/* scale thermal strain */
	if (constraint == SolidMaterialT::kPlaneStrain)
		return 1.0 + (fC12/(fC12 + fC11));
	else
		 return 1.0;
}
