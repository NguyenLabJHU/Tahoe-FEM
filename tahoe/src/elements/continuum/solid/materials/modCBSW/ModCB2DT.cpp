/* $Id: ModCB2DT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/31/1997)                                          */

#include "ModCB2DT.h"

#include <math.h>
#include <iostream.h>

#include "Constants.h"

#include "fstreamT.h"
#include "ModCBSolverT.h"
#include "dMatrixT.h"

/* material parameters */
const int knsd = 2;

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* plane codes - for crystal axes rotated wrt global axes*/
const int	kDC2Dnatural 	= 0;
const int 	kDC2D110		= 1;
const int	kDC2D111		= 2;

/* constructor */
ModCB2DT::ModCB2DT(ifstreamT& in, const ElasticT& element, bool equilibrate):
	NL_E_Mat2DT(in, element, kPlaneStrain),
	fModCBSolver(NULL),
	fCij3D(dSymMatrixT::NumValues(3)),
	fXsi(3), fStretch3D(3),
	fStretch2D(2), fStress3D(3)
{
	/* lattice transformation */
	dMatrixT Q;
	in >> fPlaneCode;
	switch (fPlaneCode)
	{
		case kDC2Dnatural:
			
			//no Q to construct
			break;

		case kDC2D110:
		{
			Q.Allocate(3);
			Q = 0.0;
			
			double cos45 = 0.5*sqrt2;
			
			/* transform global xy-plane into [110] */			
			Q(0,0) = 1.0;
			Q(1,1) = Q(2,2) = cos45;
			Q(1,2) =-cos45;
			Q(2,1) = cos45;

			break;
		}
		case kDC2D111:
		{
			Q.Allocate(3);
			Q = 0.0;
			
			/* transform global xy-plane into [111] */			
			double rt2b2 = sqrt2/2.0;
			double rt3b3 = sqrt3/3.0;
			double rt6b6 = (sqrt2*sqrt3)/6.0;
			double rt23  = sqrt2/sqrt3;
			
			Q(0,0) =-rt2b2;
			Q(0,1) =-rt6b6;
			Q(0,2) = rt3b3;
			
			Q(1,0) = rt2b2;
			Q(1,1) =-rt6b6;
			Q(1,2) = rt3b3;
			
			Q(2,0) = 0.0;
			Q(2,1) = rt23;
			Q(2,2) = rt3b3;

			break;
		}
		default:

			cout << "\nModCB2DT::ModCB2DT: unknown plane code:" << fPlaneCode;
			cout << endl;
			throw eBadInputValue;
	}

	fModCBSolver = new ModCBSolverT(Q, fThermal, in, equilibrate);
	if (!fModCBSolver) throw(eOutOfMemory);
}

/* destructor */
ModCB2DT::~ModCB2DT(void)
{
	delete fModCBSolver;
}

/* I/O functions */
void ModCB2DT::Print(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::Print(out);

	/* potential data */
	fModCBSolver->Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

void ModCB2DT::PrintName(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::PrintName(out);
	
	const char* planes[] = {"natural", "110", "111"};
	out << "    Modified CB <" << planes[fPlaneCode] << "> Plane Strain\n";

	/* potential name */
	fModCBSolver->PrintName(out);
}

void ModCB2DT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* convert 2D Green strain to 3D stretch */
	StrainToStretch(E, fStretch3D);

	/* compute moduli */
	fModCBSolver->SetModuli(fStretch3D, fXsi, fCij3D);
	
	/* convert to 2D */
	moduli.Rank4ReduceFrom3D(fCij3D);
}

/*
* Returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector.
*/
void ModCB2DT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* convert 2D Green strain to 3D stretch */
	StrainToStretch(E, fStretch3D);

	/* compute stress */
	fModCBSolver->SetStress(fStretch3D, fXsi, fStress3D);
	
	/* 3D tensor2 -> 2D vector */
	PK2[0] = fStress3D(0,0);
	PK2[1] = fStress3D(1,1);
	PK2[2] = fStress3D(0,1);	
}

/* strain energy density for the specified strain */
double ModCB2DT::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* convert 2D Green strain to 3D stretch */
	StrainToStretch(E, fStretch3D);

	return( fModCBSolver->StrainEnergyDensity(fStretch3D,fXsi) );
}

/*************************************************************************
* Private
*************************************************************************/

/*
* Compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain)
*/
void ModCB2DT::StrainToStretch(const dSymMatrixT& strain2D,
	dMatrixT& stretch3D)
{
	stretch3D = 0.0;
	
	stretch3D(0,0) = 1.0 + 2.0*strain2D[0];
	stretch3D(1,1) = 1.0 + 2.0*strain2D[1];
	stretch3D(2,2) = 1.0;
	stretch3D(0,1) = stretch3D(1,0) = 2.0*strain2D[2];	
}
