/* $Id: SSKStV2D.cpp,v 1.5 2002-11-14 17:06:06 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV2D.h"
#include "StringT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"phi", "J2_dev", "p"};

/* constructor */
SSKStV2D::SSKStV2D(ifstreamT& in, const SSMatSupportT& support):
	SSKStV(in, support),
	Material2DT(in)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* print parameters */
void SSKStV2D::Print(ostream& out) const
{
	/* inherited */
	SSKStV::Print(out);
	Material2DT::Print(out);
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables */
int SSKStV2D::NumOutputVariables(void) const { return kNumOutput; }
void SSKStV2D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
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
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void SSKStV2D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli2D(modulus, fConstraintOption);
	modulus *= fThickness;
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
		double factor = IsotropicT::DilatationFactor2D(fConstraintOption);
		thermal_strain.PlusIdentity(factor*fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}
