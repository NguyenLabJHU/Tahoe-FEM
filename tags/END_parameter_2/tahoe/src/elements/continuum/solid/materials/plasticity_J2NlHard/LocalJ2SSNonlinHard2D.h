/* $Id: LocalJ2SSNonlinHard2D.h,v 1.2.32.1 2004-03-03 16:15:02 paklein Exp $ */
#ifndef _LOCAL_J2_SS_NONLIN_HARD_2D_H_
#define _LOCAL_J2_SS_NONLIN_HARD_2D_H_

/* base classes */
#include "LocalJ2SSNonlinHard.h"

namespace Tahoe {

class LocalJ2SSNonlinHard2D : public LocalJ2SSNonlinHard
{
public:

	/* constructor */
  	LocalJ2SSNonlinHard2D(ifstreamT& in, const SSMatSupportT& support);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);

	/* print parameters */
	virtual void PrintName(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

private:

	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT    fModulus2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _LOCAL_J2_SS_NONLIN_HARD_2D_H_ */
