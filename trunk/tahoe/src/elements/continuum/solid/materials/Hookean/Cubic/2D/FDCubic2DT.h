/* $Id: FDCubic2DT.h,v 1.8 2005-01-13 00:11:24 paklein Exp $ */
/* created: paklein (06/11/1997) */
#ifndef _FD_CUBIC_2D_T_H_
#define _FD_CUBIC_2D_T_H_

/* base classes */
#include "FDCubicT.h"
#include "Anisotropic2DT.h"

namespace Tahoe {

class FDCubic2DT: public FDCubicT, public Anisotropic2DT
{
public:

	/** constructor */
	FDCubic2DT(void);

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

	/** set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
};

} /* namespace Tahoe */

#endif /* _FD_CUBIC_2D_T_H_ */
