/* $Id: QuadLogOgden2DT.h,v 1.1.2.1 2001-06-22 14:18:07 paklein Exp $ */
/* created: paklein (02/18/2001)                                          */
/* plane strain QuadLog with Ogden principal stretch formulation          */

#ifndef _QUAD_LOG_ISO_2D_T_H_
#define _QUAD_LOG_ISO_2D_T_H_

/* base classes */
#include "OgdenIsotropicT.h"
#include "Material2DT.h"

class QuadLogOgden2DT: public OgdenIsotropicT, public Material2DT
{
public:

	/* constructor */
	QuadLogOgden2DT(ifstreamT& in, const FiniteStrainT& element);

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

private:

	/* log strain */
	dArrayT flogE;		
};

#endif /* _QUAD_LOG_ISO_2D_T_H_ */
