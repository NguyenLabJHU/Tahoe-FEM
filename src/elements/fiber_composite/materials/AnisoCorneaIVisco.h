/* $Id: AnisoCorneaIVisco.h,v 1.3 2006-11-12 18:28:36 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _AnisoCorneaIVisco_
#define _AnisoCorneaIVisco_ 

/* base class */
#include "AnisoCorneaVisco.h"

#ifdef VIB_MATERIAL

namespace Tahoe {
/*forward declarations*/
class CirclePointsT;
class C1FunctionT;

class AnisoCorneaIVisco: public AnisoCorneaVisco
{
   public:
  
	/* constructor/destructor */
	AnisoCorneaIVisco(void);
		
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

protected:


	/*compute the algorithmic moduli dSNEQ/dCv deltaCv/deltadC in local fiber coord sys*/
	virtual void ComputeCalg (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v,  dMatrixT& Calg, const int process_index); 
	
	/*local newton loop for viscous stretch tensor*/ 
	virtual void Compute_Cv(const dSymMatrixT& C, const dSymMatrixT& Cv_last, dSymMatrixT& Cv, const int process_index);
				

protected:
	
	/*inverse of the viscosity function*/
	double fieta;
	dMatrixT fMod3;
	
};
	
}
#endif /* _AnisoCorneaIVisco_ */
#endif /*VIB_MATERIAL*/
