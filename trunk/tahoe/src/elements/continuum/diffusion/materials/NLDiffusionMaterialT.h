/* $Id: NLDiffusionMaterialT.h,v 1.4 2005-01-07 02:16:03 paklein Exp $ */
#ifndef _NL_DIFFUSION_MATERIALT_H_
#define _NL_DIFFUSION_MATERIALT_H_

/* base class */
#include "DiffusionMaterialT.h"

namespace Tahoe {

/* forward declarations */
class C1FunctionT;

/** interface for nonlinear materials for diffusion. Materials may
 * have temperature dependent conductivity and thermal capacity and
 * must return the associated derivatives */
class NLDiffusionMaterialT: public DiffusionMaterialT
{
public:

	/** constructor */
	NLDiffusionMaterialT(void);

	/** destructor */
	~NLDiffusionMaterialT(void);

	/** \name parameters at the current field point */
	/*@{*/
	/** conductivity */
	virtual const dMatrixT& k_ij(void);

	/** heat flux */
	virtual const dArrayT& q_i(void);
	
	/** change in heat flux with temperature */
	virtual const dArrayT& dq_i_dT(void);

	/** change in conductivity with temperature */
	virtual const dMatrixT& dk_ij(void);

	/** specific heat */
	virtual double SpecificHeat(void) const;

	/** change in capacity with temperature */
	virtual double dCapacity_dT(void) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists
	 * \param sub_lists description of subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** return the description of the given inline subordinate parameter list.
	 * Method will be called for each subordinate defined as inline by ParameterInterfaceT::SubNames
	 * or defined recursively by ParameterInterfaceT::DefineInlineSub. 
	 * \param sub name of the inlined subordinate list
	 * \param order defines whether list is a sequence or choice.
	 * \param sub_listss description of contents of this sub list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

private:

	/** \name temperature dependence */
	/*@{*/
	/** variation of conductivity with temperature */
	C1FunctionT* fConductivityScaleFunction;

	/** variation of specific heat with temperature */
	C1FunctionT* fCpScaleFunction;
	/*@}*/

	/** temperature varying conductivity return value */
	dMatrixT fScaledConductivity;
};

} /* namespace Tahoe */

#endif /* _NL_DIFFUSION_MATERIALT_H_ */
