/* $Id: ModCB3DT.cpp,v 1.8 2004-06-17 07:41:03 paklein Exp $ */
/* created: paklein (10/14/1998) */
#include "ModCB3DT.h"

#include <math.h>
#include <iostream.h>

#include "toolboxConstants.h"

#include "ifstreamT.h"
#include "ModCBSolverT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* material parameters */
const int kNSD  = 3;
const int kNDOF = 3;

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* plane codes - for crystal axes rotated wrt global axes*/
const int	kDCnatural 	= 0;
const int 	kDC110		= 1;
const int	kDC111		= 2;

/* constructor */
ModCB3DT::ModCB3DT(ifstreamT& in, const FSMatSupportT& support, bool equilibrate):
	NL_E_MatT(in, support),
	fModCBSolver(NULL),
	fXsi(kNDOF),
	fC(kNSD),
	fPK2(kNSD)
{
	/* lattice transformation */
	dMatrixT Q;
	in >> fOrientationCode;
	switch (fOrientationCode)
	{
		case kDCnatural:
			
			//no Q to construct
			break;

		case kDC110:
		{
			Q.Dimension(3);
			Q = 0.0;
			
			double cos45 = 0.5*sqrt2;
			
			/* transform global xy-plane into [110] */			
			Q(0,0) = 1.0;
			Q(1,1) = Q(2,2) = cos45;
			Q(1,2) =-cos45;
			Q(2,1) = cos45;

			break;
		}
		case kDC111:
		{
			Q.Dimension(3);
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

			cout << "\nModCB3DT::ModCB3DT: unknown orientation code:" << fOrientationCode;
			cout << endl;
			throw ExceptionT::kBadInputValue ;
	}


	fModCBSolver = new ModCBSolverT(Q, fThermal, in, equilibrate);
	if (!fModCBSolver) throw ExceptionT::kOutOfMemory;
}

/* destructor */
ModCB3DT::~ModCB3DT(void) { delete fModCBSolver; }

/* I/O functions */
void ModCB3DT::Print(ostream& out) const
{
	/* inherited */
	NL_E_MatT::Print(out);

	/* potential data */
	fModCBSolver->Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

void ModCB3DT::PrintName(ostream& out) const
{
	/* inherited */
	NL_E_MatT::PrintName(out);
	
	const char* orients[] = {"natural", "110", "111"};
	out << "    Modified CB <" << orients[fOrientationCode] << ">\n";

	/* potential name */
	fModCBSolver->PrintName(out);
}

void ModCB3DT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute moduli */
	fModCBSolver->SetModuli(fC, fXsi, moduli);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void ModCB3DT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute stress */
	fModCBSolver->SetStress(fC, fXsi, fPK2);
	
	/* shape change */
	PK2.FromMatrix(fPK2);
}

/* returns the strain energy density for the specified strain */
double ModCB3DT::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	return fModCBSolver->StrainEnergyDensity(fC, fXsi);
}

/*************************************************************************
* Private
*************************************************************************/

/* compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain */
void ModCB3DT::StrainToStretch(const dSymMatrixT& E, dMatrixT& C)
{
	/* shape change */
	E.ToMatrix(C);

	/* convert */
	C *= 2.0;
	C.PlusIdentity();
}
