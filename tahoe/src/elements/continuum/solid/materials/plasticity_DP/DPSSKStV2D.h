/* $Id: DPSSKStV2D.h,v 1.11 2004-09-10 22:39:27 paklein Exp $ */
/* created: myip (06/01/1999) */
#ifndef _DP_SS_KSTV_2D_H_
#define _DP_SS_KSTV_2D_H_

/* base class */
#include "DPSSKStV.h"

namespace Tahoe {

class DPSSKStV2D: public DPSSKStV
{
  public:

	/** constructor */
	DPSSKStV2D(void);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
				const ElementCardT& element, int ip);

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
  	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
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
#endif /* _DP_SS_KSTV_2D_H_ */
