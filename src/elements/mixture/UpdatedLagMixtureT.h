/* $Id: UpdatedLagMixtureT.h,v 1.1 2004-11-05 22:53:49 paklein Exp $ */
#ifndef _UPDATED_LAG_MIXTURE_T_H_
#define _UPDATED_LAG_MIXTURE_T_H_

/* base class */
#include "UpdatedLagrangianT.h"

namespace Tahoe {

/** update Lagrangian, finite strain solid mixture */
class UpdatedLagMixtureT: public UpdatedLagrangianT
{
public:

	/** constructor */
	UpdatedLagMixtureT(const ElementSupportT& support);

	/** resolve the species name into the index */
	int SpeciesIndex(const StringT& name) const;

	/** project the Cauchy stress for the given species to the nodes */
	void ProjectPartialStress(int i);

protected:

};

} /* namespace Tahoe */

#endif /* _UPDATED_LAG_MIXTURE_T_H_ */
