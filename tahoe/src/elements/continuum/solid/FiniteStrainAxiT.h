/* $Id: FiniteStrainAxiT.h,v 1.2 2004-02-03 08:24:57 paklein Exp $ */
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
	FiniteStrainAxiT(const ElementSupportT& support, const FieldT& field);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

  protected:

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param nsd number of spatial dimensions
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(int nsd, int size);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

  protected:
 
	/** 2D tensor workspace */
	dMatrixT fMat2D;       

	/** current coords with local ordering */
	LocalArrayT fLocCurrCoords;	

	/** \name radius to integration points computed during FiniteStrainAxiT::SetGlobalShape */
	/*@{*/
	dArrayT fRadius_X;
	dArrayT fRadius_x;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FINITE_STRAIN_AXI_T_H_ */
