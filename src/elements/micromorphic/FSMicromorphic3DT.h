/* $Id: FSMicromorphic3DT.h,v 1.21 2009-07-28 22:38:47 isbuga Exp $ */
//DEVELOPMENT
#ifndef _FS_MICROMORPHIC_3D_T_H_
#define _FS_MICROMORPHIC_3D_T_H_

/* base classes */
#include "ElementBaseT.h"
#include "StringT.h"
#include "Traction_CardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"

#include "ModelManagerT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"

#include "iAutoArrayT.h"
#include "ScheduleT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"

namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class Traction_CardT;
class StringT;

/** FSMicromorphic3DT: This class contains a coupled finite deformation
 * micromorphic (displacement and micro-displacement gradient dofs) finite
 * element implementation in 3D.  The model description can be
 * found in Regueiro ASCE JEM 135:178-191 for pressure-sensitive elastoplasticity.
 **/

class FSMicromorphic3DT: public ElementBaseT
{

public:
/* material parameters */
	enum fMaterial_T 	{
	    kMu,
	    kLambda,
	    kRho_0,
	    kg,
	    kg1,
	    kg2,
	    kg3,
	    //add to this list
	    kNUM_FMATERIAL_TERMS	};

//	enum fIntegrate_T 	{
//	    kBeta,
//	    kGamma,
//	    kNUM_FINTEGRATE_TERMS	};

	/** constructor */
	FSMicromorphic3DT(	const ElementSupportT& support );

	/** destructor */
	~FSMicromorphic3DT(void);

	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunctionDispl(void) const;
	const ShapeFunctionT& ShapeFunctionMicro(void) const;

	/** echo input */
	void Echo_Input_Data (void);

	/** return true if the element contributes to the solution of the
	 * given group. ElementBaseT::InGroup returns true if group is the
	 * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
	virtual bool InGroup(int group) const;

	/* initialize/finalize time increment */
	/*@{*/
	virtual void InitStep(void);
	virtual void CloseStep(void);
	//virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state

	/** element level reconfiguration for the current time increment */
	//virtual GlobalT::RelaxCodeT RelaxSystem(void);
	/*@}*/

	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations( AutoArrayT<const iArray2DT*>& eq_d,
				AutoArrayT<const RaggedArray2DT<int>*>& eq_phi);

	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);

	/** \name writing output */
	/*@{*/
	/** register element for output */
	virtual void RegisterOutput(void);

	/** write element output */
	virtual void WriteOutput(void);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	/*@}*/

	/** return geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, iArrayT& num_facet_nodes) const;

	/** return the geometry code */
	virtual GeometryT::CodeT GeometryCode(void) const;

	/*set active elements*/
	//virtual void SetStatus(const ArrayT<ElementCardT::StatusT>& status);

	/** initial condition/restart functions (per time sequence) */
	virtual void InitialCondition(void);

	/** mass types */
	enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2, /**< diagonally lumped mass */
             kAutomaticMass = 3  /**< select the mass type base on the time integration scheme */};
	MassTypeT static int2MassTypeT(int i);

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
				     SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void);

	void Select_Equations ( const int &iBalLinMom, const int &iBalFirstMomMom );

private:

	/** \name solution methods.
	 * Both of these drivers assemble the LHS as well as the residual.
	 */
	/*@{*/
	/** driver for staggered solution */
	void RHSDriver_staggered(void);

	/** driver for monolithic solution */
	void RHSDriver_monolithic(void);
	/*@}*/

protected:

	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;

private:

	/** Gradients and other matrices */
	dMatrixT fgrad_u, fgrad_u_n;
	dMatrixT fgrad_Phi, fgrad_Phi_n;

	dMatrixT fShapeDispl, fShapeDisplGrad, fShapeDisplGrad_t;
	dMatrixT fShapeDisplGrad_t_Transpose, fShapeDisplGradGrad, fShapeDisplGrad_temp;
	dMatrixT fShapeMicro,fShapeMicroGrad;

	dMatrixT fDefGrad, fDefGradInv, fDefGradInvMatrix;

	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode_displ, fGeometryCodeSurf_displ,
	    fGeometryCode_micro, fGeometryCodeSurf_micro;
	int fGeometryCode_displ_int, fGeometryCodeSurf_displ_int;

	/** number of integration points */
	int	fNumIP_displ, fNumIPSurf_displ, fNumIP_micro, fNumIPSurf_micro;
	int knum_d_state, knum_i_state, knumstress, knumstrain;
	int num_sidesets;

	/*@}*/

	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT u;		//solid displacement
	LocalArrayT u_n; 	//solid displacement from time t_n
	LocalArrayT u_dot; 	//solid velocity
	LocalArrayT u_dot_n; 	//solid velocity from time t_n
	LocalArrayT u_dotdot; 	//solid acceleration
	LocalArrayT u_dotdot_n; 	//solid acceleration from time t_n
	LocalArrayT del_u;	//displacement increment including Newton-R update, i.e. del_u = u_{n+1}^{k+1} - u_n
	LocalArrayT Phi;	//micro-displacement-gradient dof
	LocalArrayT Phi_dot;	//micro-displacement-gradient first time derivative
	LocalArrayT Phi_dot_n;	//micro-displacement-gradient first time derivative from time t_n
	LocalArrayT Phi_dotdot;	//micro-displacement-gradient second derivative
	LocalArrayT Phi_dotdot_n;	//micro-displacement-gradient second derivative from time t_n
	LocalArrayT Phi_n;	//micro-displacement-gradient from time t_n
	LocalArrayT del_Phi;	//micro-displacement-gradient increment
	dArrayT		del_u_vec;  	// vector form
	dArrayT		del_Phi_vec;	// vector form
	dArrayT		u_vec;  	// solid displacement in vector form
	dArrayT		u_dot_vec;  	// solid velocity in vector form
	dArrayT		u_dotdot_vec;  	// solid acceleration in vector form
	dArrayT		Phi_vec;	// micro-displacement-gradient in vector form
	dArrayT		Phi_dot_vec;	// first derivative of micro-displacement-gradient in vector form
	dArrayT		Phi_dotdot_vec;	// second derivative of micro-displacement-gradient in vector form

	/*@}*/

	// problem size definitions
	int n_en_displ, n_en_displ_x_n_sd, n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd,n_en_micro_x_n_sd_x_n_sd, n_en_micro_x_n_sd;
	int n_el, n_sd, n_sd_surf, n_en_surf;
	int n_en_micro, ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro, ndof_per_nd_micro_x_n_sd;

	int step_number;
	int iConstitutiveModelType;

	//name of output vector
	StringT output;

	dArrayT fForces_at_Node;
	bool bStep_Complete;

 	double time, fRho, fRho_0;

 	void Get_Fd_ext ( dArrayT &fFd_ext );

	//-- Material Parameters
	dArrayT fMaterial_Params;

	//-- Newmark Time Integration Parameters
	dArrayT fIntegration_Params;

	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the
	  * reference coordinates */
	ShapeFunctionT* fShapes_displ;
	ShapeFunctionT* fShapes_micro;

	/** reference coordinates */
	LocalArrayT fInitCoords_displ, fInitCoords_micro;
	/** current coordinates */
	LocalArrayT fCurrCoords_displ, fCurrCoords_micro;
	/*@}*/

	/* Data Storage */
	ElementMatrixT fKdd, fKdphi;
	ElementMatrixT fKphid, fKphiphi;
	dArrayT 	fFd_int;
	dArrayT 	fFd_ext;
	dArrayT		fFphi_int;
	dArrayT		fFphi_ext;

	dArrayT		fGrad_disp_vector;
	dArrayT 	fDefGradInv_vector;
	dArrayT 	fKirchhoff_vector;
	dArrayT 	fSecond_Piola_vector;
	dArrayT		fChi_temp_vector;
	dArrayT		fFd_int_N1_vector;
	dArrayT		fTemp_vector_ndof_se;
	dArrayT		fTemp_vector_nen_micro;
	dArrayT		fTemp_vector_9x1;
	dArrayT		fPi_temp_transpose_vector;
	dArrayT		fGrad_1_J_vector;
	dArrayT	        fTemp_nsd_vector;
	dArrayT	        fFd_int_smallstrain_vector;
	dArrayT	        fGravity_vector;
	dArrayT	        fFd_int_G4_vector;
	dArrayT         fTemp_six_values;
	dArrayT         fGradv_vector;
	dArrayT         fgradv_vector;

	dMatrixT	fDeformation_Gradient;
	dMatrixT	fDefGradT_9x9_matrix;
	dMatrixT	fRight_Cauchy_Green_tensor;
	dMatrixT	fRight_Cauchy_Green_tensor_Inverse;
	dMatrixT	fLeft_Cauchy_Green_tensor;
	dMatrixT	fLeft_Cauchy_Green_tensor_Inverse;
	dMatrixT	fDeformation_Gradient_Inverse;
	dMatrixT	fDeformation_Gradient_Transpose;
	dMatrixT	fDeformation_Gradient_Inverse_Transpose;
	dMatrixT	fDefGradInv_Grad_grad;
	dMatrixT	fDefGradInv_Grad_grad_Transpose;
	dMatrixT	fIdentity_matrix;
        dMatrixT	fSecond_Piola_tensor;
        dMatrixT	fTemp_matrix_nsd_x_nsd;
        dMatrixT	fTemp_matrix_nsd_x_1;
        dMatrixT	fTemp_matrix_ndof_se_x_ndof_se;
        dMatrixT	fTemp_matrix1_ndof_se_x_ndof_se;
        dMatrixT	fTemp_matrix_ndof_se_x_nen_micro;
        dMatrixT	fTemp_matrix1_nen_press_x_ndof_se;
        dMatrixT	fTemp_matrix_nsd_x_ndof_se;
        dMatrixT	fTemp_matrix_nsd_x_nen_micro;
        dMatrixT	fKirchhoff_tensor;
        dMatrixT	fIota_temp_matrix;
        dMatrixT	fVarpi_temp_matrix;


        dMatrixT	fIm_temp_matrix;
        dMatrixT	fHbar_temp_matrix;
        dMatrixT	fEll_temp_matrix;
        dMatrixT	fPi_temp_row_matrix;
        dMatrixT	fK_dd_G3_1_matrix;
        dMatrixT	fK_dd_G3_2_matrix;
        dMatrixT	fK_dd_G3_3_matrix;
        dMatrixT	fK_dd_G3_4_matrix;
        dMatrixT	fK_dd_G4_matrix;
        dMatrixT	fI_ij_column_matrix;
        dMatrixT	fShapeMicro_row_matrix;
        dMatrixT	fWp_temp_matrix;
        dMatrixT	fChi_temp_column_matrix;
        dMatrixT	fc_matrix;
        dMatrixT	fC_matrix;
        dMatrixT	fIm_Prim_temp_matrix;
        dMatrixT	fB_matrix;
        dMatrixT	fD_matrix;
        dMatrixT	fK_dd_BTDB_matrix;
	dMatrixT 	fDefGradInv_column_matrix;
	dMatrixT 	fDefGradInv_column_matrix_Transpose;
	dMatrixT	u_dotdot_column_matrix;
	dMatrixT	fXi_temp_matrix;
	dMatrixT	fVarsigma_temp_matrix;
	dMatrixT	fI_ijkl_matrix;
	dMatrixT	u_dot_column_matrix;
	dMatrixT	u_dot_column_matrix_Transpose;
	dMatrixT	fGravity_column_matrix;
	dMatrixT	fAleph_temp_matrix;
	dMatrixT	micro_dot_column_matrix;
	dMatrixT	fImath_temp_matrix;
	dMatrixT	fPf_0_matrix;
	//////////////////////////////////////////////////////////
	/////DEFINITIONS FOR MICROMORPHIC MATRICES////////////////
	//////////////////////////////////////////////////////////
	double KrDelta[3][3];

	//Varitional Matrices coming from the Balance of linear Momentum
	dMatrixT Tsigma_1;
	dMatrixT fG1_1;
	double SigN[3][3];
	double Fn[3][3];
	double FnInv[3][3];
	double Finv[3][3];
	double Chi[3][3];
	double ChiInv[3][3];
	double ChiN[3][3];
	dMatrixT fIota_temp_1_matrix;
	dMatrixT Finv_w; // to create Iota_1 which is different than Iota because sequence of the components in wk,l
	dMatrixT Tsigma_2;
	dMatrixT fG1_2;//not being calculated yet
    dMatrixT Tsigma_3;
    dMatrixT fG1_3;// not being calculated yet
    dMatrixT TFn_1;// to be multiplied by (lamda+Tau)
    dMatrixT fG1_4;
    dMatrixT TFn_2;// to be multiplied by (Mu+sigma)
    dMatrixT fG1_5;
    dMatrixT TFn_3;// to be multiplied by (Mu+sigma)
    dMatrixT fG1_6;
    dMatrixT TChi_1;// to be multiplied by eta
    dMatrixT fG1_7;
    dMatrixT TFn_4;// to be multiplied by eta
    dMatrixT fG1_8;
    dMatrixT TChi_2;// to be multiplied by kappa
    dMatrixT fG1_9;
    dMatrixT TFn_5;// to be multiplied by kappa
    dMatrixT fG1_10;
    dMatrixT TChi_3;// to be multiplied by nu
    dMatrixT fG1_11;
    dMatrixT TFn_6;// to be multiplied by nu
    dMatrixT fG1_12;
    dMatrixT SigCurr;
    dMatrixT fG1_13;
 // Variational Matrices coming from the Balance of First Moment of Momentum
	double mn[3][3][3];
	double Mnplus1[3][3][3];
	double GammaN[3][3][3];
	double CCof[3][3][3][3][3][3];
	double GRAD_ChiN[3][3][3];// GRADIENT  in reference configuration!
	double GRAD_Chi[3][3][3];// GRADIENT in reference configuration!

	dMatrixT Mm_1;
    dMatrixT Mm_2;
    dMatrixT Mm_3;
    dMatrixT Mm_4;
    dMatrixT Mm_5;
    dMatrixT Mm_6;
    dMatrixT Mm_7;
    dMatrixT Mm_8;
    dMatrixT Mm_9;
    dMatrixT Mm_10;
    dMatrixT Mm_11;
    dMatrixT Mm_12;
    dMatrixT Mm_13;
    dMatrixT Mm_14;
    dMatrixT Ru_1;//u
    dMatrixT Ru_2;
    dMatrixT Ru_3;
    dMatrixT RChi_1;
    dMatrixT Ru_4;
    dMatrixT RChi_2;
    dMatrixT Ru_5;
    dMatrixT Rs_sigma;
    dMatrixT R_Capital_Gamma_Chi;
    dMatrixT CapitalGamma;
    dMatrixT sn_sigman;
    dMatrixT GRAD_NCHI;
    dMatrixT Finv_eta;

    dMatrixT fH1_1;
    dMatrixT fH1_2;
    dMatrixT fH1_3;
    dMatrixT fH1_4;
    dMatrixT fH1_5;
    dMatrixT fH1_6;
    dMatrixT fH1_7;
    dMatrixT fH1_8;
    dMatrixT fH1_9;
    dMatrixT fH1_10;
    dMatrixT fH1_11;
    dMatrixT fH1_12;
    dMatrixT fH1_13;
    dMatrixT fH1_14;

    dMatrixT fH2_1;
    dMatrixT fH2_2;
    dMatrixT fH2_3;
    dMatrixT fH2_4;
    dMatrixT fH2_5;
    dMatrixT fH2_6;
    dMatrixT fH2_7;
    dMatrixT fH2_8;

    dMatrixT fH3_1;


	//////////////////////////////////////////////////////////
	/////DEFINITIONS FINISH HERE FOR MICROMORPHIC MATRICES////
	//////////////////////////////////////////////////////////

        /* to store fEulerian_strain_tensor_current_IP */
        dMatrixT	fEulerian_strain_tensor_current_IP;
        /* to store fEulerian_strain_IPs for each of the 27 IPs of each element */
        dArray2DT	fEulerian_strain_IPs;
        /* to store fCauchy_stress_tensor_current_IP */
        dMatrixT	fCauchy_stress_tensor_current_IP;
        /* to store fCauchy_stress_IPs for each of the 27 IPs of each element */
        dArray2DT	fCauchy_stress_IPs;
        dArray2DT	fState_variables_IPs;
        dArray2DT	fEulerian_strain_Elements_IPs;
        dArray2DT	fCauchy_stress_Elements_IPs;
        dArray2DT	fState_variables_Elements_IPs;


	/** the solid displacement field */
	const FieldT* fDispl;

	/** the micro-displacement-gradient field */
	const FieldT* fMicro;

	/** \name state variable storage *
	 * State variables are handled ABAQUS-style. For every iteration, the state
	 * variables from the previous increment are passed to the element, which
	 * updates the values in place. Each row in the array is the state variable
	 * storage for all integration points for an element */
	/*@{*/
	dArray2DT fdState_new;
	dArray2DT fdState;

	iArray2DT fiState_new;
	iArray2DT fiState;
	/*@}*/

	/** \name connectivities */
	/*@{*/
	ArrayT<const iArray2DT*> fConnectivities_displ;
	ArrayT<const iArray2DT*> fConnectivities_micro;
	ArrayT<iArray2DT> fConnectivities_reduced;
	/*@}*/

	/** \name equations */
	/*@{*/
	ArrayT<iArray2DT> fEqnos_displ;
	ArrayT<iArray2DT> fEqnos_micro;
	/*@}*/

	/** \name element cards */
	/*@{*/
	AutoArrayT<ElementCardT> fElementCards_displ;
	AutoArrayT<ElementCardT> fElementCards_micro;
	/*@}*/

	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;

	/** integration point stresses. Calculated and stored during
	 * FSMicromorphic2DT::RHSDriver */
	dArray2DT fIPVariable;
	/*@}*/

	/** \name prescribed plastic gradient side set ID */
	/*@{*/
	ArrayT<StringT> fSideSetID;

	/** prescribed micro-displacement-gradient weight over the side set;
	    the direction is defined by {n1,n2,n3} ?? */
	//ArrayT<iArray2DT> fMicroWght;

	/** for each side set, the global nodes on the faces in the set */
	ArrayT<iArray2DT> fMicroFaces;

	/** equation numbers for the nodes on each face */
	ArrayT<iArray2DT> fMicroFaceEqnos;

	/** side set elements */
	ArrayT<iArrayT> fSideSetElements;

	/** side set faces */
	ArrayT<iArrayT> fSideSetFaces;
	/*@}*/

	/** write output for debugging */
	/*@{*/
	/** output file stream */
	ofstreamT fs_micromorph3D_out;

	/** line output formating variables */
	int outputPrecision, outputFileWidth;
	/*@}*/

	void Form_solid_shape_functions(const double* &shapes_displ_X);
	void Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp);
	void Form_Gradient_t_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp);
	void Form_micro_shape_functions(const double* &shapes_micro_X);
	void Form_deformation_gradient_tensor(void);
	void Form_Grad_grad_transformation_matrix(void);
	void Form_fDefGradT_9x9_matrix(void);
	void Form_deformation_gradient_inv_vector(void);
	void Form_kirchhoff_stress_vector(void);
	void Form_Varpi_temp_matrix(void);
	void Form_Im_temp_matrix(void);
	void Form_Hbar_temp_matrix(void);
	void Form_Ell_temp_matrix(void);
	void Form_C_matrix(const double& J_Prim);
	void Form_c_matrix(void);
	void Form_Im_Prim_temp_matrix(void);
	void Form_D_matrix(void);
	void Form_B_matrix(void);
	void Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values);
	void Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);
	void Form_gradv_vector(void);
	void Form_Xi_temp_matrix(void);
	void Form_Varsigma_temp_matrix(void);
	void Form_I_ijkl_matrix(void);
	void Compute_norm_of_array(double& norm,const LocalArrayT& B);
	//////////////////////////////////////////////////////////
	/////FUNCTIONS  FOR MICROMORPHIC MATRICES////////////////
	//////////////////////////////////////////////////////////
	//Forming the Matrices coming from the Balance of Linear Momentum
	void Form_double_Finv_from_Deformation_tensor_inverse(void);
	void Form_Finv_w_matrix(void);
	void Form_KroneckerDelta_matrix(void);
	void Form_Tsigma_1_matrix(void);
	void Form_fG1_1_matrix(void);//not defined yet
	void Form_Tsigma_2_matrix(void);
    void Form_fG1_2_matrix(void);
    void Form_Tsigma_3_matrix(void);
    void Form_fG1_3_matrix(void);
    void Form_TFn_1_matrix(void);
    void Form_fG1_4_matrix(void);
    void Form_TFn_2_matrix(void);
    void Form_fG1_5_matrix(void);
    void Form_TFn_3_matrix(void);
    void Form_fG1_6_matrix(void);
    void Form_TChi_1_matrix(void);
    void Form_fG1_7_matrix(void);
    void Form_TFn_4_matrix(void);
    void Form_fG1_8_matrix(void);
    void Form_TChi_2_matrix(void);
    void Form_fG1_9_matrix(void);
    void Form_TFn_5_matrix(void);
    void Form_fG1_10_matrix(void);
    void Form_TChi_3_matrix(void);
    void Form_fG1_11_matrix(void);
    void Form_TFn_6_matrix(void);
    void Form_fG1_12_matrix(void);
    void Form_SigCurr_vector_Cauchy_Stss_matrix(void);
    void Form_fG1_13_matrix(void);
    //Forming the Matrices coming from the Balance of First Moment of Momentum
	void Form_Gradient_of_micro_shape_functions(const dMatrixT &fShapeMicroGrad);
    void Form_Mm_1_matrix(void);// needs to be multiplied by "-" and J
    void Form_Mm_2_matrix(void);// needs to be multiplied by J
    void Form_Mm_3_matrix(void);// needs to be multiplied by J
    void Form_Mm_4_matrix(void);// needs to be multiplied by J
    void Form_Mm_5_matrix(void);// needs to be multiplied by J
    void Form_Mm_6_matrix(void);// needs to be multiplied by J
    void Form_Mm_7_matrix(void);
    void Form_Mm_8_matrix(void);
    void Form_Mm_9_matrix(void);
    void Form_Mm_10_matrix(void);
    void Form_Mm_11_matrix(void);
    void Form_Mm_12_matrix(void);
    void Form_Mm_13_matrix(void);
    void Form_Mm_14_matrix(void);//something should be changed in the loop due to div(du)!!! changed done ok!
    void Form_Ru_1_matrix(void);
    void Form_Ru_2_matrix(void);
    void Form_Ru_3_matrix(void);
    void Form_RChi_1_matrix(void);// needs to be multiplied by Kappa and J
    void Form_Ru_4_matrix(void);// needs to be multiplied by Kappa and J
    void Form_RChi_2_matrix(void);//needs to be multiplied by Nu and J
    void Form_Ru_5_matrix(void);// needs to be multiplied by Nu and J
    void Form_Rs_sigma_matrix(void);
    void Form_R_Capital_Gamma_Chi_matrix(void);// DO NOT multiply with J !!!
    void Form_Finv_eta_matrix(void);

	//////////////////////////////////////////////////////////
	/////FUNCTIONS FINISH HERE FOR MICROMORPHIC MATRICES////
	//////////////////////////////////////////////////////////

protected:

	/** extract natural boundary condition information */
	void TakeNaturalBC(const ParameterListT& list);

	/** apply traction boundary conditions to displacement equations */
	void ApplyTractionBC(void);

  	/** update traction BC data for displacement equations */
	void SetTractionBC(void);

	/* traction data */
	ArrayT<Traction_CardT> fTractionList;
	int fTractionBCSet;

	/** \name arrays with local ordering */
	/*@{*/
	LocalArrayT fLocInitCoords;   /**< initial coords with local ordering */
	LocalArrayT fLocDisp;	      /**< solid displacements with local ordering  */
	/*@}*/

	/** \name work space */
	/*@{*/
	dArrayT fNEEvec; /**< work space vector: [element DOF] */
	dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
	/*@}*/

};


inline const ShapeFunctionT& FSMicromorphic3DT::ShapeFunctionDispl(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes_displ)
	    ExceptionT::GeneralFail("FSMicromorphic3DT::ShapeFunctionDispl", "no displ shape functions");
#endif
	return *fShapes_displ;
}

inline const ShapeFunctionT& FSMicromorphic3DT::ShapeFunctionMicro(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes_micro)
	    ExceptionT::GeneralFail("FSMicromorphic3DT::ShapeFunctionMicro", "no micro shape functions");
#endif
	return *fShapes_micro;
}

/* return the geometry code */
inline GeometryT::CodeT FSMicromorphic3DT::GeometryCode(void) const
{ return fGeometryCode_displ; }


} // namespace Tahoe
#endif /* _FS_MICROMORPHIC_3D_T_H_ */



