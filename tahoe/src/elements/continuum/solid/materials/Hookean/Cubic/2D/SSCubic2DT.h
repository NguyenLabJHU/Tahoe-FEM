/* $Id: SSCubic2DT.h,v 1.1.1.1.2.1 2001-06-06 16:22:03 paklein Exp $ */
/* created: paklein (06/11/97)                                            */

#ifndef _SS_CUBIC_2D_T_H_
#define _SS_CUBIC_2D_T_H_

/* base classes */
#include "SSCubicT.h"
#include "Material2DT.h"

class SSCubic2DT: public SSCubicT, public Material2DT
{
public:

	/* constructor */
	SSCubic2DT(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);
};

#endif /* _SS_CUBIC_2D_T_H_ */
