/* $Id: MixtureSpeciesT.h,v 1.5 2005-01-14 00:20:30 paklein Exp $ */
#ifndef _MIXTURE_SPECIES_T_H_
#define _MIXTURE_SPECIES_T_H_

/* base class */
#include "NLDiffusionElementT.h"

namespace Tahoe {

class UpdatedLagMixtureT;

/** class to handle transport with component of mixture */
class MixtureSpeciesT: public NLDiffusionElementT
{
public:
	
	/** constructor */
	MixtureSpeciesT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** method used to compute stress gradient */
	enum GradientOptionT {
		kGlobalProjection,
		kElementProjection
	};

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** compute the mass flux and flux velocities */
	void ComputeMassFlux(bool compute_dmass_flux);

	/** compute the flux velocities and their variation with concentration */
//	void ComputeDMassFlux(void);

	/** compute the divergence tensor field given the values at the integration points */
	void ComputeDivergence(const dMatrixT& ip_grad_transform, const ArrayT<dMatrixT>& tensor_ip,
		dArrayT& div) const;

protected:

	/** method used to compute stress gradient */
	GradientOptionT fGradientOption;

	/** background solid */
	UpdatedLagMixtureT* fUpdatedLagMixture;

	/** background species */
	MixtureSpeciesT* fBackgroundSpecies;
	
	/** index of the species within the mixture */
	int fIndex;

	/** flux velocities */
	dArray2DT fFluxVelocity;

	/** mass flux */
	dArray2DT fMassFlux;

	/** variation in mass flux with concentration */
	dArray2DT fDMassFlux;

	/** nodal stresses */
	dArray2DT fP_avg;
	
	/** \name element-by-element stress projection */
	/*@{*/
	ArrayT<dMatrixT> fP_ip;
	ArrayT<dMatrixT> fdP_ip;
	ArrayT<dMatrixT> fip_gradient;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MIXTURE_SPECIES_T_H_ */
