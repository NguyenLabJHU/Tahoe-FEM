/* $Id: J2SSKStV2D.h,v 1.5 2004-07-15 08:28:54 paklein Exp $ */
/* created: paklein (06/18/1997) */
#ifndef _J2_SS_KSTV_2D_H_
#define _J2_SS_KSTV_2D_H_

/* base classes */
#include "J2SSKStV.h"

namespace Tahoe {

class J2SSKStV2D: public J2SSKStV
{
public:

	/** constructor */
	J2SSKStV2D(void);

	/** returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int nip, int ip);
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* return values */
	dSymMatrixT	fStress2D;
	dMatrixT	fModulus2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _J2_SS_KSTV_2D_H_ */
