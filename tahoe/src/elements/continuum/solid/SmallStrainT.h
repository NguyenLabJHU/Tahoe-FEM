/* $Id: SmallStrainT.h,v 1.9.2.2 2002-11-13 08:40:42 paklein Exp $ */
#ifndef _SMALL_STRAIN_T_H_
#define _SMALL_STRAIN_T_H_

/* base class */
#include "ElasticT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** Interface for linear strain deformation and field gradients */
class SmallStrainT: public ElasticT
{
  public:
      
	/** constructor */
	SmallStrainT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~SmallStrainT(void);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/

	/** total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;
	/*@}*/

	/** TEMPORARY. Need this extra call here to set the source for the iteration number
	 * in SmallStrainT::fSSMatSupport. The solvers are not constructed when the material
	 * support is initialized */
	virtual void InitialCondition(void);

  protected:

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is responsible for
	 * for freeing the pointer. */
	virtual MaterialListT* NewMaterialList(int size);

	/** construct list of materials from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** initialize local field arrays. Allocate B-bar workspace if needed. */
	virtual void SetLocalArrays(void);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kstrain = 0,
	                kstrain_last = 1};

	/** compute mean shape function gradient, Hughes (4.5.23) */
	void SetMeanGradient(dArray2DT& mean_gradient) const;

  private:
    
	/** offset to material needs */
	int fNeedsOffset; //NOTE - better to have this or a separate array?
  
  	/** \name return values */
	/*@{*/
  	ArrayT<dSymMatrixT> fStrain_List;
  	ArrayT<dSymMatrixT> fStrain_last_List;
	/*@}*/
  	
  	/** \name work space */
  	/*@{*/
  	dMatrixT fGradU;
  	dArrayT fLocDispTranspose; /**< used for B-bar method */
	dArray2DT fMeanGradient;   /**< store mean shape function gradient */
  	/*@}*/

  	/** the material support used to construct materials lists. This pointer
  	 * is only set the first time SmallStrainT::NewMaterialList is called. */
	SSMatSupportT* fSSMatSupport;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
	}
#endif

	return fStrain_last_List[ip];
}

} // namespace Tahoe 
#endif /* _SMALLSTRAIN_T_H_ */
