/* $Id: SmallStrainAxiT.h,v 1.2 2004-07-15 08:26:27 paklein Exp $ */
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
	SmallStrainAxiT(const ElementSupportT& support);

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

	/** return true */
	virtual bool Axisymmetric(void) const { return true; };

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

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
