/* $Id: SimoIso2D.h,v 1.2 2001-07-03 01:35:14 paklein Exp $ */
/* created: paklein (03/04/1997)                                          */
/* (2D <-> 3D) translator for the SimoIso3D.                              */

#ifndef _SIMO_ISO_2D_H_
#define _SIMO_ISO_2D_H_

/* base classes */
#include "SimoIso3D.h"
#include "Material2DT.h"

class SimoIso2D: public SimoIso3D, public Material2DT
{
public:

	/* constructor */
	SimoIso2D(ifstreamT& in, const FiniteStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stresses */
	virtual const dSymMatrixT& s_ij(void);

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT    fModulus2D;
	
	/* workspace */
	dSymMatrixT fb_2D;		 	 	
};

#endif /* _SIMO_ISO_2D_H_ */
