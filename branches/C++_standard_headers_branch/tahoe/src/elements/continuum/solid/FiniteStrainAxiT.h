/* $Id: FiniteStrainAxiT.h,v 1.5 2004-07-15 08:26:27 paklein Exp $ */
#ifndef _FINITE_STRAIN_AXI_T_H_
#define _FINITE_STRAIN_AXI_T_H_

/* base class */
#include "FiniteStrainT.h"

namespace Tahoe {

/** finite strain, axisymmetric solid */
class FiniteStrainAxiT: public FiniteStrainT
{
  public:
      
	/** constructor */
	FiniteStrainAxiT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

  protected:

	/** indicate elements are axisymmetric */
	virtual bool Axisymmetric(void) const { return true; };

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

  protected:
 
	/** 2D tensor workspace */
	dMatrixT fMat2D;       

	/** current coords with local ordering */
	LocalArrayT fLocCurrCoords;	

	/** \name radius to integration points computed during FiniteStrainAxiT::SetGlobalShape */
	/*@{*/
	/** integration point radius in undeformed configuration */
	dArrayT fRadius_X;

	/** integration point radius in current configuration */
	dArrayT fRadius_x;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FINITE_STRAIN_AXI_T_H_ */
