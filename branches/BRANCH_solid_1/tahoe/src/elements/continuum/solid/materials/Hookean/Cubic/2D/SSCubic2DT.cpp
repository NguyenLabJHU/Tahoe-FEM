/* $Id: SSCubic2DT.cpp,v 1.1.1.1.2.1 2001-06-06 16:22:03 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#include "SSCubic2DT.h"
#include "ThermalDilatationT.h"

/* constructor */
SSCubic2DT::SSCubic2DT(ifstreamT& in, const ElasticT& element):
	SSCubicT(in, element),
	Material2DT(in)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* print parameters */
void SSCubic2DT::Print(ostream& out) const
{
	/* inherited */
	SSCubicT::Print(out);
	Material2DT::Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set (material) tangent modulus */
void SSCubic2DT::SetModulus(dMatrixT& modulus)
{
	CubicT::ComputeModuli2D(modulus, fConstraintOption);
	modulus *= fThickness;
}

/*************************************************************************
* Private
*************************************************************************/

/* set the internal thermal strain */
bool SSCubic2DT::SetThermalStrain(dSymMatrixT& thermal_strain)
{
	thermal_strain = 0.0;
	if (fThermal->IsActive())
	{
		double factor = CubicT::DilatationFactor2D(fConstraintOption);
		thermal_strain.PlusIdentity(factor*fThermal->PercentElongation());
		return true;
	}
	else
		return false;
}
