/* $Id: SSEnhLocCraigT.h,v 1.2 2004-09-24 00:52:12 cfoster Exp $ */
#ifndef _SMALL_STRAIN_ENH_LOC_CF_T_H_
#define _SMALL_STRAIN_ENH_LOC_CF_T_H_

/* base class */
#include "SmallStrainT.h"
#include "SolidElementT.h"
#include "BandT.h"

#include "HookeanMatT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** Interface for linear strain deformation and field gradients */
 class SSEnhLocCraigT: public SmallStrainT //, public HookeanMatT
{
  public:
      
	/** constructor */
	SSEnhLocCraigT(const ElementSupportT& support);

	/** destructor */
	//~SSEnhLocCraigT(void);

	/** \name total strain */
	/*@{*/
	//	const dSymMatrixT& LinearStrain(void) const;
	//	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/

	/** total strain from the end of the previous time step */
	/*@{*/
	//	const dSymMatrixT& LinearStrain_last(void) const;
	//	const dSymMatrixT& LinearStrain_last(int ip) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list. */
	/*virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	  SubListT& sub_lists) const; */

	/** return the description of the given inline subordinate parameter list */
//virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, 
				ParameterListT& mat_params) const;

protected:

	/** strain-displacement options. */
	enum StrainOptionT {kStandardB = 0, /**< standard strain-displacement matrix */
	                  kMeanDilBbar = 1  /**< mean dilatation for near incompressibility */ };

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kstrain = 0,
	                kstrain_last = 1};

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

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

  private:

	/** compute mean shape function gradient, Hughes (4.5.23) */
	//void SetMeanGradient(dArray2DT& mean_gradient) const;

  protected:
    
	bool isLocalizedTemp;
	bool isLocalized;
	BandT *fBand;
	double fH_Delta;
	bool fNoBandDilation;
	double fLocalizedFrictionCoeff;
	double fJumpIncrement;
	dMatrixT fInitialModulus;


	/** driver for calculating output values */
	/* Used to check localization - is there a more appropriate fn? */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
				   const iArrayT& e_codes, dArray2DT& e_values);

//move to surface mat model?
	dSymMatrixT FormdGdSigma(int ndof);
	dSymMatrixT FormGradActiveTensorFlowDir(int ndof);
	bool IsElementLocalized();
	void ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs);
	dArrayT Centroid();

#if 0
	/** offset to material needs */
	int fNeedsOffset; //NOTE - better to have this or a separate array?
  
	/** form of B matrix */
	StrainOptionT fStrainDispOpt;
  
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
  	 * is only set the first time SSEnhLocCraigT::NewMaterialList is called. */
	SSMatSupportT* fSSMatSupport;

#endif

};

#if 0

/* inlines */

inline const dSymMatrixT& SSEnhLocCraigT::LinearStrain(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
		ExceptionT::GeneralFail("SSEnhLocCraigT::LinearStrain", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_List[CurrIP()];
}

inline const dSymMatrixT& SSEnhLocCraigT::LinearStrain(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
		ExceptionT::GeneralFail("SSEnhLocCraigT::LinearStrain", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_List[ip];
}

inline const dSymMatrixT& SSEnhLocCraigT::LinearStrain_last(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
		ExceptionT::GeneralFail("SSEnhLocCraigT::LinearStrain_last", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_last_List[CurrIP()];
}

inline const dSymMatrixT& SSEnhLocCraigT::LinearStrain_last(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
		ExceptionT::GeneralFail("SSEnhLocCraigT::LinearStrain_last", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_last_List[ip];
}

#endif

} // namespace Tahoe 

#endif /* _SMALLSTRAIN_ENHLOC_CF_T_H_ */
