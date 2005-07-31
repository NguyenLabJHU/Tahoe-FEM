/* $Id: SW2BodyT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/20/1997)                                          */

#include "SW2BodyT.h"
#include <math.h>
#include "SWDataT.h"
#include "ThermalDilatationT.h"

/* constructor */
SW2BodyT::SW2BodyT(const dArrayT& lengths, const ThermalDilatationT* thermal,
	const SWDataT& SW):
	TwoBodyT(lengths, thermal),
	fSW(SW)
{

}

/* set free dof - triggers recomputation */
void SW2BodyT::Set(void)
{	
	double* pPhi   = fPhi.Pointer();
	double* pdPhi  = fdPhi.Pointer();
	double* pddPhi = fddPhi.Pointer();
	
	double* pl = fLengths.Pointer();
	
	/* expansion factor */
	double a = fSW.fa;
	if (fThermal) a *= (1.0 + fThermal->PercentElongation());
	
	for (int i = 0; i < fLengths.Length(); i++)
	{
		*pPhi++   = U2body(*pl,a);
		*pdPhi++  = DU2body(*pl,a);
		*pddPhi++ = DDU2body(*pl,a);
	
		pl++;
	}
}

/**********************************************************************
* Private
**********************************************************************/

	/* 2 body potential and derivatives */
double SW2BodyT::U2body(double r, double a) const
{
	return( fSW.fA*fSW.feps*(-1 + pow(a,4)*fSW.fB/pow(r,4))*
	         ( (r < fSW.frcut*a) ?
	           exp(fSW.fdelta/(-fSW.frcut + r/a)) : 0.0 ) );
}

double SW2BodyT::DU2body(double r, double a) const
{
	return( -4*pow(a,4)*fSW.fA*fSW.fB*fSW.feps*
	          ( (r < fSW.frcut*a) ?
	          	exp(fSW.fdelta/(-fSW.frcut + r/a)) : 0.0 )/pow(r,5) +
	fSW.fA*fSW.feps*(-1 + pow(a,4)*fSW.fB/pow(r,4))*
( (r < fSW.frcut*a) ?
	-(a*fSW.fdelta*exp(fSW.fdelta/(-fSW.frcut + r/a))/
				pow(-(a*fSW.frcut) + r,2)) : 0.0) );
}

double SW2BodyT::DDU2body(double r, double a) const
{
	return( 20*pow(a,4)*fSW.fA*fSW.fB*fSW.feps*
	         ( (r < fSW.frcut*a) ?
	         	exp(fSW.fdelta/(-fSW.frcut + r/a)) : 0.0)/
	         	pow(r,6) - 8*pow(a,4)*fSW.fA*fSW.fB*fSW.feps*
( (r < fSW.frcut*a) ?
	-(a*fSW.fdelta*exp(fSW.fdelta/(-fSW.frcut + r/a))/
				pow(-(a*fSW.frcut) + r,2)) : 0.0 )/pow(r,5) +
fSW.fA*fSW.feps*(-1 + pow(a,4)*fSW.fB/pow(r,4))*
( (r < fSW.frcut*a) ?
	a*fSW.fdelta*(a*fSW.fdelta - 2*a*fSW.frcut +
	2*r)*exp(fSW.fdelta/(-fSW.frcut + r/a))/
				pow(-(a*fSW.frcut) + r,4): 0.0 ) );
}
