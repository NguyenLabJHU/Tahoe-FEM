/* $Id: SimoQ1P0Axi_inv.h,v 1.2 2004-06-26 18:39:04 paklein Exp $ */
#ifndef _SIMO_Q1_P0_AXI_INV_H_
#define _SIMO_Q1_P0_AXI_INV_H_

/* base classes */
#include "UpdatedLagrangianAxiT.h"

namespace Tahoe {

/** finite strain, mixed element formulation. Using (dilation)^-1 instead
 * of the dilation. */
class SimoQ1P0Axi_inv: public UpdatedLagrangianAxiT
{
public:

	/** constructor */
	SimoQ1P0Axi_inv(const ElementSupportT& support, const FieldT& field);

	/** data initialization */
	virtual void Initialize(void);

	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;

protected:

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/** read materials data. For two dimensional problems, checks
	 * that all materials are plane strain. */
	virtual void ReadMaterialData(ifstreamT& in);	
	
private:

	/** compute mean shape function gradient, V (reference volume), and
	 * the inverse dilation */
	void SetMeanGradient(dArray2DT& mean_gradient, double& V, double& Gamma) const;

	/** special mixed index term in the tangent. Needed to compute
	 * the term in the tangent resulting from
	 * \f$ \nabla \mathbf{u} \textrm{:} \left( \nabla \boldsymbol{\eta} \right)^T \f$.
	 */
	void bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const;
	
protected:

	/** \name element volume */
	/*@{*/
	/** reference element volume */
	dArrayT fElementVolume;

	/** inverse dilation */
	dArrayT fGamma;

	/** inverse dilation from the previous time increment */
	dArrayT fGamma_last;
	/*@}*/
	
	/** element pressure. Calculated during SimoQ1P0Axi_inv::FormKd. */
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
#endif /* _SIMO_Q1_P0_AXI_INV_H_ */
