/* $Id: SimoQ1P0Axi.h,v 1.3 2004-07-15 08:26:27 paklein Exp $ */
#ifndef _SIMO_Q1_P0_AXI_H_
#define _SIMO_Q1_P0_AXI_H_

/* base classes */
#include "UpdatedLagrangianAxiT.h"

namespace Tahoe {

/** axisymmetric implementation of SimoQ1P0 */
class SimoQ1P0Axi: public UpdatedLagrangianAxiT
{
public:

	/** constructor */
	SimoQ1P0Axi(const ElementSupportT& support);

	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

private:

	/** compute mean shape function gradient, H (reference volume), and
	 * current element volume, equation (2.20) */
	void SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const;
	
	/** special mixed index term in the tangent. Needed to compute
	 * the term in the tangent resulting from
	 * \f$ \nabla \mathbf{u} \textrm{:} \left( \nabla \boldsymbol{\eta} \right)^T \f$.
	 */
	void bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const;
	
protected:

	/** \name element volume */
	/*@{*/
	/** deformed element volume */
	dArrayT fElementVolume;

	/** deformed element volume from the last time step */
	dArrayT fElementVolume_last;
	/*@}*/
	
	/** element pressure. Calculated during SimoQ1P0Axi::FormKd. */
	dArrayT fPressure;

	/** determinant of the deformation gradient for the current element */
	dArrayT fJacobian;

	/** \name work space */
	/*@{*/
	dArray2DT fMeanGradient; /**< mean gradient over element */
	dMatrixT fF_tmp; /**< F workspace */
	dMatrixT fNEEmat; /**< dimension of stiffness matrix */
	dMatrixT fdiff_b;
	dMatrixT fb_bar;
	dMatrixT fb_sig;
	/*@}*/

	/** debugging flags */
	bool fOutputInit;
	
	//TEMP - cell tracking
	int fOutputCell;
};

} /* namespace Tahoe */

#endif /* _SIMO_Q1_P0_AXI_H_ */
