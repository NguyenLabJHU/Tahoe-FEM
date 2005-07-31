/* $Id: SSKStV2D.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#include "SSKStV2D.h"
#include "StringT.h"
#include "ThermalDilatationT.h"

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"phi", "J2_dev", "p"};

/* constructor */
SSKStV2D::SSKStV2D(ifstreamT& in, const ElasticT& element):
	SSHookeanMatT(in, element),
	KStV2D(in, fModulus, fDensity)
{

}

/* print parameters */
void SSKStV2D::Print(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::Print(out);
	KStV2D::Print(out);
}

/* print name */
void SSKStV2D::PrintName(ostream& out) const
{
	/* inherited */
	SSHookeanMatT::PrintName(out);
	KStV2D::PrintName(out);
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables */
int SSKStV2D::NumOutputVariables(void) const { return kNumOutput; }
void SSKStV2D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Allocate(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void SSKStV2D::ComputeOutput(dArrayT& output)
{
	/* phi */
	output[0] = StrainEnergyDensity();

	/* compute Cauchy stress */
	const dSymMatrixT& cauchy_2D = s_ij();
	
	/* elastic constants */
	double nu = Poisson();
	
	/* 3D stress tensor */
	double a[6];
	dSymMatrixT cauchy_3D(3,a);
	cauchy_3D.ExpandFrom2D(cauchy_2D);
		
	cauchy_3D(2,2) = (fConstraintOption == kPlaneStress) ? 0.0:
						nu*(cauchy_3D(0,0) + cauchy_3D(0,0));

	/* pressure */
	output[2] = cauchy_3D.Trace()/3.0;
		
	/* deviator J2 */
	cauchy_3D.Deviatoric();
	output[1] = cauchy_3D.Invariant2();
}

/*************************************************************************
* Private
*************************************************************************/

/* set the internal thermal strain */
bool SSKStV2D::SetThermalStrain(dSymMatrixT& thermal_strain)
{
	thermal_strain = 0.0;
	if (fThermal->IsActive())
	{
		double factor = DilatationFactor();
		thermal_strain.PlusIdentity(factor*fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}
