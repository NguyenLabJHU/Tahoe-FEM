/* $Id: MixtureSpeciesT.h,v 1.2 2004-11-07 17:08:48 paklein Exp $ */
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

	/** form the residual force vector */
	virtual void RHSDriver(void);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

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
};

} /* namespace Tahoe */

#endif /* _MIXTURE_SPECIES_T_H_ */
