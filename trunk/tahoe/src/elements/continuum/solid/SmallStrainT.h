/* $Id: SmallStrainT.h,v 1.3 2002-04-16 16:35:59 paklein Exp $ */

#ifndef _SMALL_STRAIN_T_H_
#define _SMALL_STRAIN_T_H_

/* base class */
#include "ElasticT.h"

/** Interface for linear strain deformation and field gradients */
class SmallStrainT: public ElasticT
{
  public:
      
	/** constructor */
	SmallStrainT(FEManagerT& fe_manager);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** total strain */
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;

	/** total strain from the end of the previous time step */
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;

	/** compute field gradients */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients from the end of the previous time step */
	void ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

  protected:

	/** construct list of materials from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** initialize local field arrays. Allocate B-bar workspace if needed. */
	virtual void SetLocalArrays(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);
  	           
  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kstrain = 0,
	                kstrain_last = 1};

  private:
    
	/** offset to material needs */
	int fNeedsOffset; //NOTE - better to have this or a separate array?
  
  	/** return values */
  	ArrayT<dSymMatrixT> fStrain_List;
  	ArrayT<dSymMatrixT> fStrain_last_List;
  	
  	/** \name work space */
  	/*@{*/
  	dMatrixT fGradU;
  	dArrayT  fLocDispTranspose; /**< used for B-bar method */
  	/*@}*/
};

/* inlines */

inline const dSymMatrixT& SmallStrainT::LinearStrain(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
	{
		cout << "\n SmallStrainT::LinearStrain: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
#endif

	return fStrain_List[CurrIP()];
}

inline const dSymMatrixT& SmallStrainT::LinearStrain(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
	{
		cout << "\n SmallStrainT::LinearStrain: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
#endif

	return fStrain_List[ip];
}

inline const dSymMatrixT& SmallStrainT::LinearStrain_last(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
	{
		cout << "\n SmallStrainT::LinearStrain_last: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
#endif

	return fStrain_last_List[CurrIP()];
}

inline const dSymMatrixT& SmallStrainT::LinearStrain_last(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
	{
		cout << "\n SmallStrainT::LinearStrain_last: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw eGeneralFail;
	}
#endif

	return fStrain_last_List[ip];
}

#endif /* _SMALLSTRAIN_T_H_ */
