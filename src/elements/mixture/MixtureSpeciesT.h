/* $Id: MixtureSpeciesT.h,v 1.3 2005-01-03 21:55:34 paklein Exp $ */
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

	/** compute the flux velocities */
	void ComputeMassFlux(void);

protected:

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
	
	/** \name finite difference of mass flux */
	/*@{*/
	dArray2DT fFluxVelocity_tmp;
	dArray2DT fDMassFlux;	
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MIXTURE_SPECIES_T_H_ */
