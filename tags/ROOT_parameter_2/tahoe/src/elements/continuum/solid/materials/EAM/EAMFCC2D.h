/* $Id: EAMFCC2D.h,v 1.7 2003-01-29 07:34:38 paklein Exp $ */
/* created: paklein (12/09/1996) */
#ifndef _EAMFCC2D_H_
#define _EAMFCC2D_H_

/* base class */
#include "NL_E_Mat2DT.h"

namespace Tahoe {

/* forward declarations */
class EAMFCC3DSym;

/** plane strain EAM material */
class EAMFCC2D: public NL_E_Mat2DT
{
public:

	/** plane codes - for crystal axes rotated wrt global axes */
	enum PlaneCodeT {kFCC001 = 0,
                     kFCC101 = 1,
                     kFCC111 = 2};

	/* constructor */
	EAMFCC2D(ifstreamT& in, const FSMatSupportT& support, PlaneCodeT plane_code);

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
	
	PlaneCodeT fPlaneCode;
	int	fEAMCode;
	
	/* EAM solver */
	EAMFCC3DSym* fEAM;
	
};

} // namespace Tahoe 
#endif /* _EAMFCC2D_H_ */
