/* $Id: LJFCC111.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (07/31/1996)                                          */

#ifndef _LJFCC111_H_
#define _LJFCC111_H_

/* base class */
#include "NL_E_RotMat2DT.h"

class LJFCC111: public NL_E_RotMat2DT
{
public:

	/* constructor */
	LJFCC111(ifstreamT& in, const ElasticT& element);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					                 					
	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

private:

	/* scaling constant */
	double fScale;
};

#endif /* _LJFCC111_H_ */
