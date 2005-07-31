/* $Id: EAMFCC2D.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (12/09/1996)                                          */
/* Plane strain EAM material                                              */

#include "EAMFCC2D.h"

#include <math.h>
#include <iostream.h>

#include "Constants.h"

#include "fstreamT.h"
#include "EAMFCC3DSym.h"
#include "dMatrixT.h"

/* material parameters */
const int knsd = 2;

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* constructor */
EAMFCC2D::EAMFCC2D(ifstreamT& in, const ElasticT& element, int planecode):
	NL_E_Mat2DT(in, element, kPlaneStrain),
	fPlaneCode(planecode),
	fEAM(NULL)
{
	/* read EAM code */
	in >> fEAMCode;
	
	switch (fPlaneCode)
	{
		case kFCC2Dnatural:
		{
			fEAM = new EAMFCC3DSym(in, fEAMCode, knsd);
			break;
		}	
		case kFCC2D110:
		{
			dMatrixT Q(3);
			
			double cos45 = 0.5*sqrt2;
			
			/* transform global xy-plane into [110] */
			Q = 0.0;
			
			Q(0,0) = 1.0;
			Q(1,1) = Q(2,2) = cos45;
			Q(1,2) =-cos45;
			Q(2,1) = cos45;
			
			fEAM = new EAMFCC3DSym(in, Q, fEAMCode, knsd);
			break;
		}	
		case kFCC2D111:
		{
			dMatrixT Q(3);
			Q = 0.0;
			
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
			
			fEAM = new EAMFCC3DSym(in, Q, fEAMCode, knsd);
			break;
		}	
		default:
		{
			cout << "\nEAMFCC2D::EAMFCC2D: unknown plane code:" << fPlaneCode;
			cout << endl;
			throw eBadInputValue;
		}
	}
	
	if (!fEAM) throw eOutOfMemory;
	
	fEAM->Initialize();	
}

/* destructor */
EAMFCC2D::~EAMFCC2D(void)
{
	delete fEAM;
}

/* I/O functions */
void EAMFCC2D::Print(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::Print(out);

	/* print EAM solver data */
	fEAM->Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

void EAMFCC2D::PrintName(ostream& out) const
{
	/* inherited */
	NL_E_Mat2DT::PrintName(out);

	const char* planes[] = {"natural", "110", "111"};

	out << "    EAM FCC 2D <" << planes[fPlaneCode] << "> Plane Strain\n";
}

/*************************************************************************
* Private
*************************************************************************/

void EAMFCC2D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* EAM solver */
	fEAM->Moduli(moduli, E);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void EAMFCC2D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* EAM solver */
	fEAM->SetStress(E, PK2);
}

/* returns the strain energy density for the specified strain */
double EAMFCC2D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* EAM solver */
	return fEAM->EnergyDensity(E);
}
