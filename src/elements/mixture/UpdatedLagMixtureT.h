/* $Id: UpdatedLagMixtureT.h,v 1.2 2004-11-07 17:09:26 paklein Exp $ */
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

	/** return the body force vector */
	void BodyForce(dArrayT& body_force) const;

	/** return the nodal accelerations over the current element */
	void Acceleration(LocalArrayT& acc);

	/** \name selecting current element externally */
	/*@{*/
	/** reset loop */
	virtual void Top(void) { UpdatedLagrangianT::Top(); };
	
	/** advance to next element. \return true if there is another element, 
	 * false otherwise */ 
	virtual bool NextElement(void) { return UpdatedLagrangianT::NextElement(); };

	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void) { UpdatedLagrangianT::SetGlobalShape(); };
	/*@}*/

protected:

};

} /* namespace Tahoe */

#endif /* _UPDATED_LAG_MIXTURE_T_H_ */
