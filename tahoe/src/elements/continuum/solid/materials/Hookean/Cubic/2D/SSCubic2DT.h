/* $Id: SSCubic2DT.h,v 1.5 2002-11-14 17:06:05 paklein Exp $ */
/* created: paklein (06/11/97) */
#ifndef _SS_CUBIC_2D_T_H_
#define _SS_CUBIC_2D_T_H_

/* base classes */
#include "SSCubicT.h"
#include "Anisotropic2DT.h"
#include "Material2DT.h"

namespace Tahoe {

class SSCubic2DT: public SSCubicT, public Anisotropic2DT, public Material2DT
{
public:

	/* constructor */
	SSCubic2DT(ifstreamT& in, const SSMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;

	/** return the pressure associated with the last call to 
	 * StructuralMaterialT::s_ij. See StructuralMaterialT::Pressure
	 * for more information. \note plane strain not implemented, but 
	 * could be using CubicT::DilatationFactor2D. */
	virtual double Pressure(void) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

} // namespace Tahoe 
#endif /* _SS_CUBIC_2D_T_H_ */
