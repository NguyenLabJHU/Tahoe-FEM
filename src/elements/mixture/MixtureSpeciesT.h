/* $Id: MixtureSpeciesT.h,v 1.1 2004-11-05 22:53:49 paklein Exp $ */
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

protected:

	/** background solid */
	UpdatedLagMixtureT* fUpdatedLagMixture;
	
	/** index of the species within the mixture */
	int fIndex;
};

} /* namespace Tahoe */

#endif /* _MIXTURE_SPECIES_T_H_ */
