/* $Id: SSCubic2DT.h,v 1.2.6.1 2002-06-27 18:03:12 cjkimme Exp $ */
/* created: paklein (06/11/97)                                            */

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
	SSCubic2DT(ifstreamT& in, const SmallStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

} // namespace Tahoe 
#endif /* _SS_CUBIC_2D_T_H_ */
