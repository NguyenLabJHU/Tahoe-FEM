/* $Id: FiniteStrainT_dev.h,v 1.2 2003-02-12 01:48:10 paklein Exp $ */
#ifndef _FINITE_STRAIN_T_H_
#define _FINITE_STRAIN_T_H_

/* base class */
#include "SolidElementT.h"

namespace Tahoe {

/* forward declarations */
class FSMatSupportT;
class FSSolidMatT;

/** Interface for linear strain deformation and field gradients */
class FiniteStrainT: public SolidElementT
{
  public:
      
	/** constructor */
	FiniteStrainT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~FiniteStrainT(void);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** TEMPORARY. Need this extra call here to set the source for the iteration number
	 * in SmallStrainT::fSSMatSupport. The solvers are not constructed when the material
	 * support is initialized */
	virtual void InitialCondition(void);

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

	/** register self for output */
	virtual void RegisterOutput(void);

	/** send output */
	virtual void WriteOutput(void);

  protected:

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is responsible for
	 * for freeing the pointer.
	 * \param size number of materials in the list */
	virtual MaterialListT* NewMaterialList(int size);

	/** construct list of materials from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

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

	/*******************************************************************/
	/*material force evaluation*/
	bool MatForceDriver(void);
	void AssembleMatForce(const dArrayT& elem_val, dArrayT& global_val, 
			      const nArrayT<int>& ndnos);
	void MatForceVolMech(FSSolidMatT* CurrMaterial, dArrayT& MatForce);
	void MatForceSurfMech(void);
	void MatForceDissip(FSSolidMatT* CurrMaterial, dArrayT& MatForce, 
			     dArrayT& statev);
	/******************************************************************/
  
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

	/*material force*/
	dArrayT fMatForce;
	dArrayT fMatForceDissip;
  
  private:
  
  	/** the material support used to construct materials lists. This pointer
  	 * is only set the first time FiniteStrainT::NewMaterialList is called. */
	FSMatSupportT* fFSMatSupport;
  
	/** offset to material needs */
	int fNeedsOffset; //NOTE - better to have this or a separate array?
	
	/** material force output ID */
	int fMatForceOutputID;
	OutputSetT fOutputSet;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
