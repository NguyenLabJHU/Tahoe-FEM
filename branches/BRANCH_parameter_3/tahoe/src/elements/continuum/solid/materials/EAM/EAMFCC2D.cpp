/* $Id: EAMFCC2D.cpp,v 1.8.46.2 2004-06-09 23:17:32 paklein Exp $ */
/* created: paklein (12/09/1996) */
#include "EAMFCC2D.h"

#include <math.h>
#include <iostream.h>

#include "toolboxConstants.h"

#include "fstreamT.h"
#include "EAMFCC3DSym.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* material parameters */
const int knsd = 2;

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* constructor */
EAMFCC2D::EAMFCC2D(ifstreamT& in, const FSMatSupportT& support, PlaneCodeT plane_code):
	ParameterInterfaceT("EAM_FCC_2D"),
	NL_E_MatT(in, support),
	fPlaneCode(plane_code),
	fEAM(NULL)
{
	/* read EAM code */
	in >> fEAMCode;
	
	switch (fPlaneCode)
	{
		case kFCC001:
		{
			fEAM = new EAMFCC3DSym(in, fEAMCode, knsd);
			break;
		}	
		case kFCC101:
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
		case kFCC111:
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
			throw ExceptionT::kBadInputValue;
		}
	}
	
	if (!fEAM) throw ExceptionT::kOutOfMemory;
	
	fEAM->Initialize();	
}

/* destructor */
EAMFCC2D::~EAMFCC2D(void)
{
	delete fEAM;
}

/* describe the parameters needed by the interface */
void EAMFCC2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
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
