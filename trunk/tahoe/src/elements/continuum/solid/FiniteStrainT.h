/* $Id: FiniteStrainT.h,v 1.11 2002-07-17 00:02:10 paklein Exp $ */

#ifndef _FINITE_STRAIN_T_H_
#define _FINITE_STRAIN_T_H_

/* base class */
#include "ElasticT.h"

namespace Tahoe {

/** Interface for linear strain deformation and field gradients */
class FiniteStrainT: public ElasticT
{
  public:
      
	/** constructor */
	FiniteStrainT(const ElementSupportT& support, const FieldT& field);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

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

  protected:

	/** return a pointer to a new material list.
	 * \param size number of materials in the list */
	virtual MaterialListT* NewMaterialList(int size) const;

	/** construct list of materials from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** calculate the damping force contribution ("-c*v"). \note arises from
	 * support for Rayleigh damping, which is on the way out */
	virtual void FormCv(double constC);

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

  	/* work space  */
  	ArrayT<dMatrixT> fF_List;      /**< deformation gradient */
  	dArrayT          fF_all;       /**< grouped memory for all deformation gradients */
  	ArrayT<dMatrixT> fF_last_List; /**< last deformation gradient */
  	dArrayT          fF_last_all;  /**< grouped memory for all last deformation gradients */

	/** Pointer to shape functions wrt current coords. This pointer must be
	 * set by sub-classes to enable calculation wrt current coordinates */
	ShapeFunctionT* fCurrShapes;
  
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
	{
		cout << "\n FiniteStrainT::DeformationGradient: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
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
	{
		cout << "\n FiniteStrainT::DeformationGradient: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
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
	{
		cout << "\n FiniteStrainT::DeformationGradient_last: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
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
	{
		cout << "\n FiniteStrainT::DeformationGradient_last: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
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
