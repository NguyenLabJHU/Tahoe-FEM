/* $Id: UpdatedLagMixtureT.h,v 1.9 2006-04-14 15:28:32 thao Exp $ */
#ifndef _UPDATED_LAG_MIXTURE_T_H_
#define _UPDATED_LAG_MIXTURE_T_H_

/* base class */
#include "UpdatedLagrangianT.h"

namespace Tahoe {

/* forward declarations */
class FSSolidMixtureT;

/** update Lagrangian, finite strain solid mixture */
class UpdatedLagMixtureT: public UpdatedLagrangianT
{
public:

	/** constructor */
	UpdatedLagMixtureT(const ElementSupportT& support);

	/** resolve the species name into the index */
	int SpeciesIndex(const StringT& name) const;

	/** density of the given species */
	double Density(int i);

	/** \name concentration flag */
	/*@{*/
	/** concentration type */
	enum ConcentrationT {
		kReference,
		kCurrent
	};
	
	/** set concentration flag 
	 * \param i species index
	 * \param conc concentration flag */
	void SetConcentration(int i, ConcentrationT conc);
	/*@}*/

	/** \name stress global projection */
	/*@{*/
	/** project the given partial first Piola-Kirchoff stress to the nodes */
	void ProjectPartialStress(int i);

	/** project the given partial first Piola-Kirchoff stress to the nodes */
	void ProjectPartialCauchy(int i);

	/** project the variation with concentration of the given partial first 
	 * Piola-Kirchoff stress to the nodes */
	void ProjectDPartialStress(int i);
	/*@}*/

	/** project the variation with concentration of the given partial cauchy stress to the nodes */
	void ProjectDPartialCauchy(int i);
	/*@}*/

	/** collect the integration point stresses and variations in stress. This method
	 * assumes all shape functions and other element quantities have been updated for
	 * the "current" element.
	 * \param i species index
	 * \param ip_stress pointer to the destination of the integration point stresses
	 *        or NULL if the integration point stresses are wanted
	 * \param ip_dstress pointer to the destination of the integration point stress
	 *        variations or NULL if the variations are wanted */
	void IP_PartialStress(int i, ArrayT<dMatrixT>* ip_stress, ArrayT<dMatrixT>* ip_dstress);

	/*returns the species specific kirchhoff stress at the integration points*/
	void IP_PartialCauchy(int i, ArrayT<dMatrixT>* ip_stress, ArrayT<dMatrixT>* ip_dstress);

	/** return the body force vector */
	void BodyForce(dArrayT& body_force) const;

	/** return the nodal accelerations over the current element */
	void Acceleration(LocalArrayT& acc) const;

	/** return the nodal velocities over the current element */
	void Velocity(LocalArrayT& vel) const;

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

	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values);

protected:

	/** \name cast the material FSSolidMixtureT */
	/*@{*/
	FSSolidMixtureT& FSSolidMixture(void);
	const FSSolidMixtureT& FSSolidMixture(void) const;
	/*@}*/
	
protected:

	/** \name work space */
	/*@{*/
	dMatrixT fF_inv;
	dMatrixT fStress;
	/*@}*/
};

inline FSSolidMixtureT& UpdatedLagMixtureT::FSSolidMixture(void) {
	const UpdatedLagMixtureT& const_this = const_cast<UpdatedLagMixtureT&>(*this);
	return const_cast<FSSolidMixtureT&>(const_this.FSSolidMixture());
}

} /* namespace Tahoe */

#endif /* _UPDATED_LAG_MIXTURE_T_H_ */
