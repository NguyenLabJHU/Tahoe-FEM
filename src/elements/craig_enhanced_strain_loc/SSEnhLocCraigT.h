/* $Id: SSEnhLocCraigT.h,v 1.14 2005-12-15 18:39:07 cfoster Exp $ */
#ifndef _SMALL_STRAIN_ENH_LOC_CF_T_H_
#define _SMALL_STRAIN_ENH_LOC_CF_T_H_

/* base class */
#include "SmallStrainT.h"
#include "SolidElementT.h"
#include "BandT.h"

#include "HookeanMatT.h"
#include "MapT.h"

#include "ofstreamT.h"

namespace Tahoe{

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

	enum BVPTypeT {kNonhomogeneous = 0,
             kHomogeneous,
             kPreFailed,
             kNumBVPTypes};

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
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

    int fFirstElementToLocalize; // -1 if loading is nonhomogeneous, element number otherwise

	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is 
	 * responsible for freeing  the pointer. 
	 * \param name list identifier \param size length of the list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);
	virtual bool LocalizationHasBegun(void) {return
						 fLocalizationHasBegun;};
	virtual void PreFailElements();
						 

  private:

	/** compute mean shape function gradient, Hughes (4.5.23) */
	//void SetMeanGradient(dArray2DT& mean_gradient) const;
	MapT<int, BandT*> fTracedElements;
	
	/*elements with one edge on band of traced elements
	 * with that coordinate of the end of band*/
	iAutoArrayT fEdgeOfBandElements;
	AutoArrayT<dArrayT> fEdgeOfBandCoords;
	static bool fLocalizationHasBegun;
	static double fDetAMin;
	static int fLeastDetEle;
	
	ofstreamT jump_out;
	
  protected:
    
	//BVPTypeT fBVPType;
	int fBVPType;
	dArrayT fPreFailedNormal;
	dArrayT fPreFailedSlipDir;

	BandT *fBand;
	double fH_delta_0;
	bool fMultiBand;
	bool fNoBandDilation;
	double fLocalizedFrictionCoeff;
	MapT<int, BandT*>* TracedElements()
	  {return &fTracedElements;}

	/** driver for calculating output values */
	/* Used to check localization - is there a more appropriate fn? */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
				   const iArrayT& e_codes, dArray2DT& e_values); 
				   
	#if 0			   
	virtual void GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
    const iArrayT& e_codes, ArrayT<StringT>& e_labels) const;
    virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
     iArrayT& counts) const;			   
	#endif	
				   
  public:
	virtual void CloseStep(void);
  protected:
	virtual void GetElement(int elementNumber);
	virtual double CalculateJumpIncrement();
	virtual bool IsBandActive();
	virtual void LoadBand(int elementNumber);
	virtual BandT* FormNewBand(dArrayT normal, dArrayT slipDir,
				   dArrayT perpSlipDir, dArrayT coords, double area);
	virtual void AddNewEdgeElements(int elementNumber);
	virtual dArrayT InterceptCoords(dArrayT& localizedEleCoord, dArrayT& nodalCoord1, dArrayT& nodalCoord2);

//move to surface mat model?
	dSymMatrixT FormdGdSigma(int ndof);
	dSymMatrixT FormGradActiveTensorFlowDir(int ndof, int ip);
	bool IsElementTraced();
	bool IsElementTraced(int elementNumber);
	bool IsElementLocalized();
	bool TraceElement();
	virtual void ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs);
	dArrayT Centroid();

};

/* These are conforming strains. Change to regular strains? */

/* inlines */
#if 0

inline const dSymMatrixT& SSEnhLocCraigT::LinearStrain(void) const
{
#if 0 
	//__option(extended_errorcheck)
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
#if 0
//__option(extended_errorcheck)
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
#if 0
//__option(extended_errorcheck)
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
#if 0 
//__option(extended_errorcheck)
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
