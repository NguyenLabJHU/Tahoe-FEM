/* $Id: SmallStrainEnhLocT.h,v 1.17 2005-04-07 00:24:26 raregue Exp $ */
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

	enum fElementLocScalars_T {
							kdetAmin,
							kdissip_max,
							kzeta,
							kgamma_delta,
							kQ_S,
							kP_S,
							kq_St,
							ksign_q_St,
							kr_S,
							kKzetazeta,
							ksecphi2,
							kh_phi,
							kDelta_zeta,
							kNUM_SCALAR_TERMS
							};
							
	enum fElementLocInternalVars_T {
							kCohesion,
							kFriction,
							kDilation,
							kNUM_ISV_TERMS
							};					
							
	enum fCohesiveSurfaceParams_T {
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

	/* check for localization */
	void CheckLocalization(int& elem, LocalArrayT& displ_elem);
	
	/* choose the normal and slipdir given normals and slipdirs from bifurcation condition */
	void ChooseNormalAndSlipDir(LocalArrayT& displ_elem, int& elem, int& nen);
	
	/* given the normal and one point, determine active nodes */
	void DetermineActiveNodesTrace(LocalArrayT& coords_elem, int& elem, int& nen);

	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);
	
	/** element level localization checks */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

private:

	/** compute mean shape function gradient, Hughes (4.5.23) */
	//void SetMeanGradient(dArray2DT& mean_gradient) const;
	/** compute mean shape function gradient, and element volume, equation (2.20) */
	void SetMeanGradient(dArray2DT& mean_gradient, double& v) const;
	
	/** write output for debugging */
	/*@{*/
	/** flag to indicate first pass, and debugging */
	//static bool fFirstPass, fDeBug, fFirstTrace;
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
	
	/** \name return values */
	/*@{*/
  	ArrayT<dSymMatrixT> fStress_List;
  	ArrayT<dSymMatrixT> fStress_last_List;
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
	dArray2DT fElementLocTangent;
	dArray2DT fElementLocSlipDir;
	dArrayT fElementLocPsi;
	dArray2DT fElementLocMuDir;
	dArray2DT fElementLocInternalVars;
	dArray2DT fElementLocGradEnh; // varies for each IP
	dArray2DT fElementLocGradEnhIP; // for each IP for one element
	dArray2DT fElementLocEdgeIntersect;
	dArray2DT fElementLocStartSurface;
	dArrayT fElementVolume;
	dArrayT fElementYieldTrial;
	
	dArray2DT fElementStress;
  	
	iArray2DT fElementLocNodesActive;
	
	iArrayT fElementLocFlag;
	
	dArray2DT fElementLastIterateDisp; // used to store last iterated displacement
	dArrayT fLocLastIterateDisp, fLocdeltaDisp; // used to calculate zeta increment
	
	// used to store last iterated K_zetad and K_dzeta
	dArray2DT fElementLocKzetad, fElementLocKdzeta;
	
	/* element centroid */
	dArray2DT fElementCentroid;
	
	double psi, psi_tmp, psi_chosen;
	
	/** from the last time step */
	dArray2DT fElementLocScalars_last;
	dArray2DT fElementLocSlipDir_last;
	dArray2DT fElementLocMuDir_last;
	dArray2DT fElementLocInternalVars_last;
	dArrayT fElementVolume_last;
	dArray2DT fElementStress_last;
	/*@}*/
	
	dArray2DT stress_IPs, stress_last_IPs;
	
	dArrayT fCohesiveSurface_Params;
							
	double c_r, c_p, alpha_c, phi_r, phi_p, alpha_phi, psi_p, alpha_psi;
	
	double detAmin, detA_tmp;
	double dissip_max, dissip_tmp;
	AutoArrayT <dArrayT> normals;
	AutoArrayT <dArrayT> slipdirs;
	AutoArrayT <dArrayT> tangents;
	AutoArrayT <double> psis;
	AutoArrayT <double> detAs;
	AutoArrayT <double> dissipations_fact;
	AutoArrayT <double> grad_displ_mns;
	AutoArrayT <dArrayT> normals_min;
	AutoArrayT <dArrayT> slipdirs_min;
	AutoArrayT <dArrayT> tangents_min;
	AutoArrayT <double> psis_min;
	AutoArrayT <double> detAs_min;
	AutoArrayT <double> dissipations_fact_min;
	AutoArrayT <double> grad_displ_mns_min;
	
	dArrayT grad_enh_IP, mu_dir, mu_dir_last;
	dArray2DT grad_enh_IPs;
		
	dArrayT normal_tmp, normal_chosen;
	dArrayT slipdir_tmp, slipdir_chosen;
	dArrayT tangent_tmp, tangent_chosen;
	dArrayT tmp_array;
	
	int loc_flag, numedges;
	dArrayT node_displ, node_coords, node_shape_deriv;
	dArrayT start_surface_vect, start_surface_vect_read, elem_centroid;
	
	int choose_normal;
	
	double fYieldTrial, r_S, K_zetazeta;
	double Q_S, Q_S_last, Q_Sn_trial, P_S, q_St, P_S_trial, q_St_trial, P_S_last, 
			q_St_last, q_St_abs, q_St_abs_trial, q_St_abs_last;
	double sign_q_St, gamma_delta, DgammadeltaDzeta, DpsiDzeta, DphiDzeta;
	double DQ_SDzeta, DQ_SnDzeta, DP_SDzeta, Dr_SDzeta;
	double h_c, h_phi, h_psi;
	double Dh_cDzeta, Dh_phiDzeta, Dh_psiDzeta;
	
	dArrayT q_isv, q_isv_last, h_q, DqDzeta, DhqDzeta, DmudirDzeta, DslipdirDzeta, DslipdirDpsi;
	
	dMatrixT fDe;
	
	dSymMatrixT F_nn, F_tn, F_mun, F_mun_last, G_enh, tmp_sym_matrix, inner_matrix;
	dSymMatrixT DF_munDzeta, DGenhDzeta;
	
	dMatrixT F_mun_nonsym, F_nn_nonsym, tmp_matrix;
	
	dSymMatrixT strain_inc, stress_inc, stress_last, stress_return;
	
	double inner, inner1, inner2, inner_abs;
	
	dSymMatrixT fStressTrial, fStressCurr, DsigDzeta, tmp_stress1, tmp_stress2;
	
	dMatrixT fK_dd, fK_dzeta_x_Kzetad;
	dArrayT fK_dzeta, fK_zetad, tmp_vec1, tmp_vec2;
	
	double zeta, zeta_last, delta_zeta, Delta_zeta;
	
	double var1, var2, var3;
	
	double cospsi, sinpsi, tanphi, tanphi_n, secphi, secphi2;

	bool fFirstTrace;

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
