/* $Id: SMRSSKStV2D.h,v 1.1 2006-07-27 13:20:08 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003) */
#ifndef _SMR_SS_KSTV_2D_H_
#define _SMR_SS_KSTV_2D_H_

/* base class */
#include "SMRSSKStV.h"

namespace Tahoe {

class SMRSSKStV2D: public SMRSSKStV
{
  public:

	/* constructor */
	SMRSSKStV2D(void);

	/* returns 3D strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
				const ElementCardT& element, int ip);
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
  	
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
  	dMatrixT	fModulus2D, fModulusPerfPlas2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _SMR_SS_KSTV_2D_H_ */