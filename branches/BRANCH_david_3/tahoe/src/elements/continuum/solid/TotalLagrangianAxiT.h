/* $Id: TotalLagrangianAxiT.h,v 1.4 2004-07-15 08:26:27 paklein Exp $ */
#ifndef _TOTAL_LAGRANGRIAN_AXI_T_H_
#define _TOTAL_LAGRANGRIAN_AXI_T_H_

/* base class */
#include "FiniteStrainAxiT.h"

namespace Tahoe {

/** total Lagrangian, finite strain element */
class TotalLagrangianAxiT: public FiniteStrainAxiT
{
public:

	/** constructors */
	TotalLagrangianAxiT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
		
protected:

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

protected:

  	/** array of shape functions */
  	dArrayT fIPShape;

	/** \name workspace */
	/*@{*/
	dMatrixT fStressMat;   /**< space for a stress 3D tensor */
	dMatrixT fStressStiff; /**< compact stress stiffness contribution */
	dMatrixT fGradNa;      /**< shape function gradients matrix */
	
	dArrayT   fTemp2;
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x;
	/*@}*/
	
	/** debugging flags */
	bool fOutputInit;
	
	//TEMP - cell tracking
	int fOutputCell;
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_AXI_T_H_ */
