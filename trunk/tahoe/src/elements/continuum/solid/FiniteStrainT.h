/* $Id: FiniteStrainT.h,v 1.2 2001-07-03 01:34:50 paklein Exp $ */

#ifndef _FINITE_STRAIN_T_H_
#define _FINITE_STRAIN_T_H_

/* base class */
#include "ElasticT.h"

/** Interface for linear strain deformation and field gradients */
class FiniteStrainT: public ElasticT
{
  public:
      
	/** constructor */
	FiniteStrainT(FEManagerT& fe_manager);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** total deformation gradient */
	const dMatrixT& DeformationGradient(void) const;
	const dMatrixT& DeformationGradient(int ip) const;

	/** total strain from the end of the previous time step */
	const dMatrixT& DeformationGradient_last(void) const;
	const dMatrixT& DeformationGradient_last(int ip) const;

	/** compute field gradients with respect to current coordinates */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients with respect to reference coordinates */
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

  protected:

	/** construct list of materials from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kF = 0,
	                kF_last = 1};

  private:
  
	/** offset to material needs */
	int fNeedsOffset; //NOTE - better to have this or a separate array?
  
  	/** return values */
  	ArrayT<dMatrixT> fF_List;
  	ArrayT<dMatrixT> fF_last_List;
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

#endif /* _FINITE_STRAIN_T_H_ */
