/* $Id: QuadLogOgden3DT.h,v 1.5 2003-01-29 07:34:43 paklein Exp $ */
/* created: paklein (02/17/2001) */
#ifndef _QUAD_LOG_ISO_3D_T_H_
#define _QUAD_LOG_ISO_3D_T_H_

/* base class */
#include "OgdenIsotropicT.h"

namespace Tahoe {

/** principal stretch version of Quad Log model */
class QuadLogOgden3DT: public OgdenIsotropicT
{
public:

	/* constructor */
	QuadLogOgden3DT(ifstreamT& in, const FSMatSupportT& support);
	
	/* print parameters */
	virtual void PrintName(ostream& out) const;

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* principal values given principal values of the stretch tensors,
	 * i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod);

protected:

	/* log strain */
	dArrayT flogE;
};

} // namespace Tahoe 
#endif /* _QUAD_LOG_ISO_3D_T_H_ */
