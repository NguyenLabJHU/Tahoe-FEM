/* created: Karma Yonten (03/04/2004)                   
   MR version modified to incorporate gradient plasticity 
   theory.
*/
#ifndef _GRAD_MR_SS_KSTV_2D_H_
#define _GRAD_MR_SS_KSTV_2D_H_

/* base class */
#include "GRAD_MRSSKStV.h"

namespace Tahoe 
{

class GRAD_MRSSKStV2D: public GRAD_MRSSKStV
{
  public:

	/* constructor */
	GRAD_MRSSKStV2D(void);
	
	/* returns 3D strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
		
	/* returns 3D Laplacian of strain (3D) */
	virtual const dSymMatrixT& LapElasticStrain(
                const dSymMatrixT& laptotalstrain, 
		const ElementCardT& element, int ip);
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
  	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);
	
	/* yield function */
	virtual const double& Yield_Function(void);
	
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
  	dMatrixT	fModulus2D, fModulusPerfPlas2D;
  	double      fYieldFunction2D; //yield function

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _GRAD_MR_SS_KSTV_2D_H_ */