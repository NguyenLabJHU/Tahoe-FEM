/* $Id: EAMFCC2D.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (12/09/1996)                                          */
/* Plane strain EAM material                                              */

#ifndef _EAMFCC2D_H_
#define _EAMFCC2D_H_

/* base class */
#include "NL_E_Mat2DT.h"

/* forward declarations */
class EAMFCC3DSym;


class EAMFCC2D: public NL_E_Mat2DT
{
public:

	/* plane codes - for crystal axes rotated wrt global axes*/
	enum PlaneCodeT {kFCC2Dnatural = 0,
                         kFCC2D110 = 1,
                         kFCC2D111 = 2};

	/* constructor */
	EAMFCC2D(ifstreamT& in, const ElasticT& element, int planecode);

	/* destructor */
	virtual ~EAMFCC2D(void);
	
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
	
protected:
	
	int fPlaneCode;
	int	fEAMCode;
	
	/* EAM solver */
	EAMFCC3DSym* fEAM;
	
};

#endif /* _EAMFCC2D_H_ */
