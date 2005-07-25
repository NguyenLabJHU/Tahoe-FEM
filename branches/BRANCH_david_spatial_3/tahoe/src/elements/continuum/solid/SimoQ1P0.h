/* $Id: SimoQ1P0.h,v 1.9.30.1 2005-07-25 02:36:55 paklein Exp $ */
#ifndef _SIMO_Q1_P0_H_
#define _SIMO_Q1_P0_H_

/* base classes */
#include "UpdatedLagrangianT.h"

namespace Tahoe {

/** finite strain, mixed element formulation.
 * The formulation is due to Simo, Taylor, and Pister, CMAME \b 51, 
 * 177-208, 1985. The basic idea behind the formulation is to
 * represent the pressure and dilatation \f$ \Theta \f$ as separate 
 * fields from the displacement. For the continuous case, the determinant 
 * of the deformation gradient 
   \f[ 
       J = \det \mathbf{F} 
         = \mathbf{1} + \frac{\partial \mathbf{u}}{\partial \mathbf{X}} 
   \f]
 * is equal to the dilatation. However, when the displacement
 * field \f$ \mathbf{u} \f$ is restricted to a finite dimensional
 * representation, it may not contain enough degrees of freedom to
 * represent nearly incompressible deformations without making the 
 * response overly stiff. Therefore, the dilatation and pressure are
 * represented as separate fields. The modified deformation gradient
 * is given by
   \f[ 
       \bar{\mathbf{F}} = \left( \frac{\Theta}{J} \right)^{1/3} \mathbf{F}.
   \f]
 * The remainder of the formulation results as a consequence. For Q1P0,
 * the (2D) bi- or (3D) trilinear displacement field (Q1) is combined
 * with piecewise constant (P0) pressure and dilatation fields. Since the
 * space for these fields is restricted to within element domains, these
 * degrees of freedom can be determined analytically at the element level
 * and substituted into the remaining element equations to results in
 * a purely displacement-based formulation.
 *
 * \note Several errors appear in the derivation
 * of the consistent tangent in the CMAME paper. Therefore,
 * the implementation of the tangent here does not match the 
 * published formulation. */
class SimoQ1P0: public UpdatedLagrangianT
{
public:

	/** constructor */
	SimoQ1P0(const ElementSupportT& support);

	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(ifstreamT& in);

	/** write restart information from stream */
	virtual void WriteRestart(ofstreamT& out) const;

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
	
	/** element pressure. Calculated during SimoQ1P0::FormKd. */
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
};

} // namespace Tahoe 
#endif /* _SIMO_Q1_P0_H_ */
