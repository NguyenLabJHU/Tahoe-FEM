/* $Id: SSCubic2DT.h,v 1.8 2005-01-13 00:11:24 paklein Exp $ */
/* created: paklein (06/11/97) */
#ifndef _SS_CUBIC_2D_T_H_
#define _SS_CUBIC_2D_T_H_

/* base classes */
#include "SSCubicT.h"
#include "Anisotropic2DT.h"

namespace Tahoe {

class SSCubic2DT: public SSCubicT, public Anisotropic2DT
{
public:

	/** constructor */
	SSCubic2DT(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. \note plane strain not implemented, but 
	 * could be using CubicT::DilatationFactor2D. */
	virtual double Pressure(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/** set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

} /* namespace Tahoe */

#endif /* _SS_CUBIC_2D_T_H_ */
