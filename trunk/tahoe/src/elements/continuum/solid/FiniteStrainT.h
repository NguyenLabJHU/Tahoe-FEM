/* $Id: FiniteStrainT.h,v 1.19 2005-03-02 02:29:14 paklein Exp $ */
#ifndef _FINITE_STRAIN_T_H_
#define _FINITE_STRAIN_T_H_

/* base class */
#include "SolidElementT.h"

namespace Tahoe {

/* forward declarations */
class FSMatSupportT;

/** Interface for linear strain deformation and field gradients */
class FiniteStrainT: public SolidElementT
{
  public:
      
	/** constructor */
	FiniteStrainT(const ElementSupportT& support);

	/** destructor */
	~FiniteStrainT(void);

	/** \name deformation gradients */
	/*@{*/
	/** total deformation gradient at the current integration point */
	const dMatrixT& DeformationGradient(void) const;

	/** total deformation gradient at the specified integration point */
	const dMatrixT& DeformationGradient(int ip) const;

	/** total strain from the end of the previous time step at the current integration point */
	const dMatrixT& DeformationGradient_last(void) const;

	/** total strain from the end of the previous time step at the specified integration point */
	const dMatrixT& DeformationGradient_last(int ip) const;
	/*@}*/

	/** \name field gradients */
	/*@{*/
	/** compute field gradients with respect to current coordinates at the current integration point */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;

	/** compute field gradients with respect to current coordinates at the specified integration point */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients with respect to reference coordinates at the current integration point */
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const;

	/** compute field gradients with respect to reference coordinates at the specified integration point */
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list. */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** return the description of the given inline subordinate parameter list */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

  protected:

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param name list identifier
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** returns true if the material requires the deformation gradient */
	bool Needs_F(int material_number) const;

	/** returns true if the material requires the deformation gradient 
	 * from the end of the last time increment */
	bool Needs_F_last(int material_number) const;

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kF = 0,
	                kF_last = 1};

  protected:

  	/** \name work space  */
  	/*@{*/
  	ArrayT<dMatrixT> fF_List;      /**< deformation gradient */
  	dArrayT          fF_all;       /**< grouped memory for all deformation gradients */
  	ArrayT<dMatrixT> fF_last_List; /**< last deformation gradient */
  	dArrayT          fF_last_all;  /**< grouped memory for all last deformation gradients */
  	/*@}*/

	/** Pointer to shape functions wrt current coords. This pointer must be
	 * set by sub-classes to enable calculation wrt current coordinates */
	ShapeFunctionT* fCurrShapes;
  
  	/** the material support used to construct materials lists. This pointer
  	 * is only set the first time FiniteStrainT::NewMaterialList is called. */
	FSMatSupportT* fFSMatSupport;

  private:
    
	/** offset to material needs */
	int fNeedsOffset; //NOTE - better to have this or a separate array?
};

/* inlines */

inline const dMatrixT& FiniteStrainT::DeformationGradient(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kF])
		ExceptionT::GeneralFail("FiniteStrainT::DeformationGradient",
			"material %d did not specify this need", mat_num+1);
#endif

	return fF_List[CurrIP()];
}

inline const dMatrixT& FiniteStrainT::DeformationGradient(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kF])
		ExceptionT::GeneralFail("FiniteStrainT::DeformationGradient",
			"material %d did not specify this need", mat_num+1);
#endif

	return fF_List[ip];
}

inline const dMatrixT& FiniteStrainT::DeformationGradient_last(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kF_last])
		ExceptionT::GeneralFail("FiniteStrainT::DeformationGradient_last",
			"material %d did not specify this need", mat_num+1);
#endif

	return fF_last_List[CurrIP()];
}

inline const dMatrixT& FiniteStrainT::DeformationGradient_last(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kF_last])
		ExceptionT::GeneralFail("FiniteStrainT::DeformationGradient_last",
			"material %d did not specify this need", mat_num+1);
#endif

	return fF_last_List[ip];
}

/* returns true if the material requires the deformation gradient */
inline bool FiniteStrainT::Needs_F(int material_number) const
{
	/* material information */
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	return needs[fNeedsOffset + kF];
}

/* returns true if the material requires the deformation gradient 
 * from the end of the last time increment */
inline bool FiniteStrainT::Needs_F_last(int material_number) const
{
	/* material information */
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	return needs[fNeedsOffset + kF_last];
}

/* compute field gradients with respect to reference coordinates at the current integration point */
inline void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const
{
	/* inherited function is wrt reference coordinates */
	IP_ComputeGradient(u, grad_u);
}

/* compute field gradients with respect to reference coordinates at the specified integration point */
inline void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{
	/* inherited function is wrt reference coordinates */
	IP_ComputeGradient(u, grad_u, ip);
}

} // namespace Tahoe 
#endif /* _FINITE_STRAIN_T_H_ */
