/* $Id: SSKStV2D.cpp,v 1.6 2004-07-15 08:27:18 paklein Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV2D.h"
#include "StringT.h"
#include "ThermalDilatationT.h"

using namespace Tahoe;

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"phi", "J2_dev", "p"};

/* constructor */
SSKStV2D::SSKStV2D(void):
	ParameterInterfaceT("small_strain_StVenant_2D")
{

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
		
	cauchy_3D(2,2) = (Constraint() == kPlaneStress) ? 0.0:
						nu*(cauchy_3D(0,0) + cauchy_3D(0,0));

	/* pressure */
	output[2] = cauchy_3D.Trace()/3.0;
		
	/* deviator J2 */
	cauchy_3D.Deviatoric();
	output[1] = cauchy_3D.Invariant2();
}

/* describe the parameters needed by the interface */
void SSKStV2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSKStV::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStress);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void SSKStV2D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli2D(modulus, Constraint());
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
		double factor = IsotropicT::DilatationFactor2D(Constraint());
		thermal_strain.PlusIdentity(factor*fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}
