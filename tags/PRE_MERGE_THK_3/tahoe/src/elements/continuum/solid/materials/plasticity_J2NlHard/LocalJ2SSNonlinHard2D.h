/* $Id: LocalJ2SSNonlinHard2D.h,v 1.2 2002-11-14 17:06:29 paklein Exp $ */
#ifndef _LOCAL_J2_SS_NONLIN_HARD_2D_H_
#define _LOCAL_J2_SS_NONLIN_HARD_2D_H_

/* base classes */
#include "LocalJ2SSNonlinHard.h"
#include "Material2DT.h"

namespace Tahoe {

class LocalJ2SSNonlinHard2D : public LocalJ2SSNonlinHard, public Material2DT
{
public:

        /* constructor */
  	LocalJ2SSNonlinHard2D(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

private:

	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT    fModulus2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _LOCAL_J2_SS_NONLIN_HARD_2D_H_ */
