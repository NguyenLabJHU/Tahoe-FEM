/* $Id: DPSSKStVLoc2D.h,v 1.7 2005-11-23 22:36:04 raregue Exp $ */
/* created: myip (06/01/1999) */
#ifndef _DP_SS_KSTV_LOC_2D_H_
#define _DP_SS_KSTV_LOC_2D_H_

/* base class */
#include "DPSSKStVLoc.h"

namespace Tahoe {

/* forward declarations */
class SSEnhLocMatSupportT;

class DPSSKStVLoc2D: public DPSSKStVLoc
{
public:

	/* constructor */
	DPSSKStVLoc2D(void);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
				const dSymMatrixT& totalstrain, 
				const ElementCardT& element, int ip);

	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& ce_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
  	
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

	// pointer to material support
	//const SSEnhLocMatSupportT* fSSEnhLocMatSupport;
  
	/* return values */
	dSymMatrixT	fStress2D;
	dMatrixT fModulus2D, fModulusElas2D, fModulusPerfPlas2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _DP_SS_KSTV_LOC_2D_H_ */
