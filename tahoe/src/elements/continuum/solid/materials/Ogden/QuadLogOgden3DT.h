/* $Id: QuadLogOgden3DT.h,v 1.2 2001-07-03 01:35:13 paklein Exp $ */
/* created: paklein (02/17/2001)                                          */
/* principal stretch version of Quad Log model                            */

#ifndef _QUAD_LOG_ISO_3D_T_H_
#define _QUAD_LOG_ISO_3D_T_H_

/* base class */
#include "OgdenIsotropicT.h"

class QuadLogOgden3DT: public OgdenIsotropicT
{
public:

	/* constructor */
	QuadLogOgden3DT(ifstreamT& in, const FiniteStrainT& element);
	
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

#endif /* _QUAD_LOG_ISO_3D_T_H_ */
