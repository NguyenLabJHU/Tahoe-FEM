/* $Id: UpdatedLagMixtureT.h,v 1.4 2005-01-14 00:20:30 paklein Exp $ */
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

	/** \name stress global projection */
	/*@{*/
	/** project the given partial first Piola-Kirchoff stress to the nodes */
	void ProjectPartialStress(int i);

	/** project the variation with concentration of the given partial first 
	 * Piola-Kirchoff stress to the nodes */
	void ProjectDPartialStress(int i);
	/*@}*/

	/** collect the integration point stresses and variations in stress.
	 * \param ip_stress pointer to the destination of the integration point stresses
	 *        or NULL if the integration point stresses are wanted
	 * \param ip_dstress pointer to the destination of the integration point stress
	 *        variations or NULL if the variations are wanted
	 * \param i */
	void IP_PartialStress(int i, ArrayT<dMatrixT>* ip_stress, ArrayT<dMatrixT>* ip_dstress);

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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name work space */
	/*@{*/
	dMatrixT fF_inv;
	dMatrixT fStress;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _UPDATED_LAG_MIXTURE_T_H_ */
