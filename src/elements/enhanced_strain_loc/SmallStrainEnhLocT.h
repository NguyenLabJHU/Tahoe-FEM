/* $Id: SmallStrainEnhLocT.h,v 1.7 2005-02-15 01:21:07 raregue Exp $ */
#ifndef _SMALL_STRAIN_ENH_LOC_T_H_
#define _SMALL_STRAIN_ENH_LOC_T_H_

/* base class */
#include "SolidElementT.h"

#include "ofstreamT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** Interface for linear strain deformation and field gradients */
class SmallStrainEnhLocT: public SolidElementT
{

public:

	enum fElementLocScalars {
							kLocFlag,
							kJumpDispl,
							kgamma_delta,
							kQ,
							kP,
							kq_St,
							kq_Sn,
							kp_S,
							kNUM_SCALAR_TERMS 
							};
							
	enum fElementLocInternalVars {
							kCohesion,
							kFriction,
							kDilation,
							kNUM_ISV_TERMS 
							};					
							
	enum fElementLocCohesiveSurfaceParams {
							kc_r,
							kc_p,
							kalpha_c,
							kphi_r,
							kphi_p,
							kalpha_phi,
							kpsi_p,
							kalpha_psi,
							kNUM_CS_TERMS 
							};												

	/** constructor */
	SmallStrainEnhLocT(const ElementSupportT& support);

	/** destructor */
	~SmallStrainEnhLocT(void);
	
	/** initialize current step */
	virtual void InitStep(void);
	
	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;

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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

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

	/* choose the normal and slipdir given normals and slipdirs from bifurcation condition */
	void ChooseNormalAndSlipDir(void);
	
	/* given the normal and one point, determine active nodes */
	void DetermineActiveNodes(void);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

private:

	/** compute mean shape function gradient, Hughes (4.5.23) */
	//void SetMeanGradient(dArray2DT& mean_gradient) const;
	/** compute mean shape function gradient, and element volume, equation (2.20) */
	void SetMeanGradient(dArray2DT& mean_gradient, double& v) const;
	
	/** write output for debugging */
	/*@{*/
	/** flag to indicate first pass, and debugging */
	static bool fFirstPass, fDeBug;
	/** output file stream */
	ofstreamT ss_enh_out;
	/** line output formating variables */
	int outputPrecision, outputFileWidth;
	/*@}*/

protected:
    
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
  	 * is only set the first time SmallStrainEnhLocT::NewMaterialList is called. */
	SSMatSupportT* fSSMatSupport;
	
/* for post-localization */
protected:

	/** \name element localization info */
	/*@{*/
	/** current time step */
	dArray2DT fElementLocScalars;
	dArray2DT fElementLocNormal;
	dArray2DT fElementLocNormal1;
	dArray2DT fElementLocNormal2;
	dArray2DT fElementLocNormal3;
	dArray2DT fElementLocTangent;
	dArray2DT fElementLocTangent1;
	dArray2DT fElementLocTangent2;
	dArray2DT fElementLocTangent3;
	dArray2DT fElementLocSlipDir;
	dArray2DT fElementLocSlipDir1;
	dArray2DT fElementLocSlipDir2;
	dArray2DT fElementLocSlipDir3;
	dArray2DT fElementLocPsiSet;
	dArray2DT fElementLocMuDir;
	dArray2DT fElementLocInternalVars;
	dArray2DT fElementLocGradEnh; // varies for each IP
	dArray2DT fElementLocGradEnhIP; // for each IP for one element
	dArray2DT fElementLocEdgeIntersect;
	dArrayT fElementVolume;
	
	double psi1, psi2, psi3;
	
	/** from the last time step */
	dArray2DT fElementLocScalars_last;
	dArray2DT fElementLocSlipDir_last;
	dArray2DT fElementLocMuDir_last;
	dArray2DT fElementLocInternalVars_last;
	dArrayT fElementVolume_last;
	/*@}*/
	
	dArrayT fCohesiveSurface_Params;
	
	AutoArrayT <dArrayT> normals;
	AutoArrayT <dArrayT> slipdirs;
	dArrayT grad_enh, mu_dir;
	dArrayT normal1, normal2, normal3, normal_chosen;
	dArrayT slipdir1, slipdir2, slipdir3, slipdir_chosen;
	dArrayT tangent1, tangent2, tangent3, tangent_chosen;
	int loc_flag, numedges;
	
	double fYieldTrial, residual_slip, K_zetazeta;
	double DgammadeltaDzeta, DpsiDzeta, DPDzeta;
	dArrayT q_isv, h_q, DqDzeta, DhqDzeta, DslipdirDzeta;
	dMatrixT F_mun, G_enh, fDe;
	dSymMatrixT F_nn;
	
	ElementMatrixT fK_dd, fK_dzeta, fK_zetad;
	
	LocalArrayT displ_u;

};

/* inlines */

inline const dSymMatrixT& SmallStrainEnhLocT::LinearStrain(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
		ExceptionT::GeneralFail("SmallStrainEnhLocT::LinearStrain", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_List[CurrIP()];
}

inline const dSymMatrixT& SmallStrainEnhLocT::LinearStrain(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
		ExceptionT::GeneralFail("SmallStrainEnhLocT::LinearStrain", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_List[ip];
}

inline const dSymMatrixT& SmallStrainEnhLocT::LinearStrain_last(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
		ExceptionT::GeneralFail("SmallStrainEnhLocT::LinearStrain_last", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_last_List[CurrIP()];
}

inline const dSymMatrixT& SmallStrainEnhLocT::LinearStrain_last(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
		ExceptionT::GeneralFail("SmallStrainEnhLocT::LinearStrain_last", 
		"material %d did not specify this need", 
			mat_num + 1);
#endif

	return fStrain_last_List[ip];
}

} // namespace Tahoe 

#endif /* _SMALLSTRAIN_ENHLOC_T_H_ */
