/* $Id: TotalLagrangianAxiT.h,v 1.2.14.2 2004-05-11 03:59:45 paklein Exp $ */
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
	TotalLagrangianAxiT(const ElementSupportT& support, const FieldT& field);
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
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_AXI_T_H_ */