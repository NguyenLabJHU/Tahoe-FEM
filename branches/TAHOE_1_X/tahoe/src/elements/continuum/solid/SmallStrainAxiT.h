/* $Id: SmallStrainAxiT.h,v 1.1 2004-01-31 07:20:48 paklein Exp $ */
#ifndef _SMALL_STRAIN_AXI_T_H_
#define _SMALL_STRAIN_AXI_T_H_

/* base class */
#include "SmallStrainT.h"

namespace Tahoe {

/** small strain, torionless, axisymmetric element */
class SmallStrainAxiT: public SmallStrainT
{
  public:
      
	/** constructor */
	SmallStrainAxiT(const ElementSupportT& support, const FieldT& field);
	SmallStrainAxiT(const ElementSupportT& support);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/*@}*/

  protected:

	/** return true */
	virtual bool Axisymmetric(void) const { return true; };

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

//TEMP - not needed with ParameterListT input
	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param nsd number of spatial dimensions
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(int nsd, int size);

	/** initialize local field arrays. Allocate B-bar workspace if needed. */
	virtual void SetLocalArrays(void);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

  private:

	/** compute mean shape function gradient, Hughes (4.5.23) */
	void SetMeanGradient(dArray2DT& mean_gradient) const;

  private:
  
  	/** space for calculating ip coordinates and displacements */
  	dArrayT fIPInterp;
  	
  	/** array of shape functions */
  	dArrayT fIPShape;
  	
  	/** 2D strain */
  	dSymMatrixT fStrain2D;

  	/** 2D-axis stress */
  	dSymMatrixT fStress2D_axi;
};

} /* namespace Tahoe */

#endif /* _SMALL_STRAIN_AXI_T_H_ */
