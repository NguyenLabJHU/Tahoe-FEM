/* $Id: StaggeredMultiScaleT.h,v 1.21 2004-07-15 08:28:27 paklein Exp $ */ 
//DEVELOPMENT
#ifndef _STAGGERED_MULTISCALE_T_H_ 
#define _STAGGERED_MULTISCALE_T_H_ 

#include "ContinuumT.h"

#define RENDER 0  // <-- Turn rendering on/off

#if RENDER
#include "Render_ManagerT.h"
#endif

/* base classes */
#include "ElementBaseT.h"

#include "StringT.h"
#include "ofstreamT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

/* base multiscale classes */
#include "FEA.h"
#include "VMS.h"
#include "FEA_FormatT.h"

namespace Tahoe {

/* forward declarations */

class ShapeFunctionT;
class Traction_CardT;	/** mass types */

/** StaggeredMultiScaleT: This class contains methods pertaining to kinematics of
 * a dual field formulation. These include deformation gradients Fe and Fp
 * and gradients such as Grad(ue) and Grad(up) as examples.
 * Sandia National Laboratory and the University of Michigan **/

class StaggeredMultiScaleT: public ElementBaseT
{
 //----- Class Methods -----------
	
 public:

	enum fMat_T 	{ 
									k__E,
									k__Pr,
									k__E1,
									k__Pr1,
									k__E2,
									k__Pr2,
									k__f,
									k__V,
									k__Y,
									k__c_zeta,
									k__l,
									k__K,
									k__H,
									k__Yield_Strain,
									k__Beta_tilde,
									k__Rho_tilde,
									k__Pi,
									k__Rho,
									k__Gamma_b,
									k__AlphaY,
									k__Density,
									kNUM_FMAT_TERMS	};		// MAT for material here, not matrix

	enum iMat_T 	{ 
									k__BS_Type,
									k__IH_Type,
									kNUM_IMAT_TERMS	};

	enum bLogic_T 	{ 
									k__Control_Eb,
									k__Diagnosis_Variables,
									k__Del_Curl_sE,
									k__BCJ_file_read,
									kNUM_LOGICAL_SWITCHES	};

	/** constructor */
	StaggeredMultiScaleT(const ElementSupportT& support, const FieldT& coarse, 
		const FieldT& fine);

	/** destructor */
	~StaggeredMultiScaleT(void);

	/** data initialization */
	virtual void Initialize(void); 

	/** echo input */
	void Echo_Input_Data (void); 

	/** return true if the element contributes to the solution of the
	 * given group. ElementBaseT::InGroup returns true if group is the
	 * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
	virtual bool InGroup(int group) const;

	/** close current time increment. Called if the integration over the
	 * current time increment was successful. */
	virtual void CloseStep(void);

	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
	
	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);

	/** register element for output */
	virtual void RegisterOutput(void);

	/** write element output */
	virtual void WriteOutput(void);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	/*@}*/

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	void Select_Equations ( const int &iCoarseScale, const int &iFineScale );

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
	
public:	
protected:
private:

	/** Data at time steps n and n+1 used by both Coarse and Fine */
	//VMS_VariableT n,np1; // <-- keep local scope in elmt loop for now 

	/** Gradients with respect to reference coodinates */
	FEA_dMatrixT fGRAD_ua, fGRAD_ua_n, fGRAD_ub, fGRAD_ub_n, fVar;
	//FEA_dScalar_ArrayT  S;

	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode;

	/** number of integration points */
	int	fNumIP;
	/*@}*/

	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT ua;     /**< fine scale displacement */
	LocalArrayT ua_n; 	/**< fine scale displacement from previous increment */
	LocalArrayT del_ua; /**< the Newton-R update i.e. del_ua = ua - ua_n (ua subcript n+1 implied) */
	LocalArrayT ub;     /**< coarse scale displacement */
	LocalArrayT DDub;     /**< coarse scale acceleration (used for body force) */
	LocalArrayT ub_n; 	/**< coarse scale displacement from previous increment */
	LocalArrayT del_ub; /**< the Newton-R update i.e. del_ub = ub - ub_n (ub subcript n+1 implied) */
	dArrayT			del_ua_vec;  	/** need in vector for i.e. { {ua1_1,ua1_2},{ua2_1,ua2_2}, ... } */
	dArrayT			del_ub_vec;		/** need in vector for i.e. { {ub1_1,ub1_2},{ub2_1,ub2_2}, ... } */
	/*@}*/

	int n_ip, n_sd, n_df, n_en, n_en_x_n_df; 
	int n_np, n_el, n_comps;
	//int step_number_last_iter;
	//bool New_Step;
	int step_number;
	int iFineScaleModelType;

	//------- Render Variables 
	
 	bool render_switch;	
	StringT render_settings_file_name;
	StringT surface_file_name;
 	double render_time;	
 	int num_tensors_to_render;	
 	int num_scalars_to_render;	
	ArrayT < StringT > Render_Tensor_Names; 
	ArrayT < StringT > Render_Scalar_Names; 
	StringT write_file_name;
	StringT node_force_file_name;
	int component_i;
	int component_j;
	int Elmt2Write;
	int ElmtIP2Write;
	bool bLog_Strain;
	int cube_top_elmt;
	int cube_top_elmt_top_local_node;
	int cube_bottom_elmt;
	int cube_bottom_elmt_bottom_local_node;
 	int render_variable_group;	
 	int render_variable_order;	
 	int render_displ;	

  ofstreamT var_plot_file;
  ofstreamT nodal_forces_file;
  ofstreamT displacements_file;
 	double x_top;  // for calc Lf	
 	double x_bot;  // for calc Lf	

	ArrayT < FEA_dMatrix_ArrayT >  Render_Tensor; // ( n_el x num variables to render (n_rv)
	ArrayT < FEA_dScalar_ArrayT >  Render_Scalar; // ( n_el x num variables to render (n_rv)

	int iDesired_Force_Node_Num, iDesired_Force_Direction;
	dArrayT fForces_at_Node;

	//-- Room for expansion w/o changing input decks
	int num_extra_integer_vars;
	int num_extra_double_vars;
	iArrayT Extra_Integer_Vars;
	dArrayT Extra_Double_Vars;

	// above will actually be ( n_el x n_rv x n_ip x n_sd x n_sd )

	bool render_data_stored;
	bool bStep_Complete;
	bool bControl_Eb, bDiagnosis_Variables, bDel_Curl_sE;
 	double time;	

	int e_set;
	FEA_dMatrixT rVariable; // Render Variable
	ContinuumT Geometry; 
	void 	Init_Render ( void );
	//void  Update_New_Step ( void );
	void 	Store_Render_Data ( double &time,LocalArrayT fInitCoords,
								 						LocalArrayT &ua,LocalArrayT &ub,FineScaleT *fEquation_2 );
	void RenderOutput ( void );
	void Get_Fext_I 	( dArrayT &fFext_1 );
	
#if RENDER
	Render_ManagerT Render_Boss;
#endif

	//-- Material Parameters 

	dArrayT fMaterial_Data;
	iArrayT iMaterial_Data;
	ArrayT <bool> bLogical_Switches;
	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the 
	 * coordinates in StaggeredMultiScaleT::fCurrCoords, which are the
	 * current coordinates */
	ShapeFunctionT* fShapes;
	
	FEA_ShapeFunctionT fFEA_Shapes;

	/** reference coordinates */
	LocalArrayT fInitCoords;     

	/** current coordinates */
	LocalArrayT fCurrCoords;
	/*@}*/

	/** the BLACK BOXS */
	/** VMF is acronym for "Variational Multi-Field". Class VMF_Virtual_WorkT is simply
	 * linearizaiton of the Virtual Work Equaiton, using a decomposiiton u = ua + ub.  This
	 * class can be used for any multi-field formulation where u is decomposed as stated, 
	 * (to include Variational Multi-Scale (VMS) formulation). */

	/* Data Storage */
	ElementMatrixT 	fKa_1, 		fKb_1;
	ElementMatrixT 	fKa_2, 		fKb_2;
	dArrayT 				fFint_1;
	dArrayT 				fFext_1;
	dArrayT					fR_2;

	/* Multi-Field Element Formulators */
	CoarseScaleT* fEquation_1;	
	FineScaleT* 	fEquation_2;

	/* Multi-Field Materials */
	VMF_MaterialT* fFineMaterial;
	VMF_MaterialT* fCoarseMaterial;

	/* Conversion methods: puts data in FEA format (very little cost in perspective) */
	FEA_FormatT Convert;

	/** the coarse scale field */
	const FieldT& fCoarse;
	
	/** the fine scale field */
	const FieldT& fFine;	

	/** equations per element for the fine scale. The coarse scale equations are
	 * in ElementBaseT::fEqnos and are handled by ElementBaseT. */
	iArray2DT fEqnos_fine;

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
	
	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;
	
	/** integration point stresses. Calculated and stored during 
	 * StaggeredMultiScaleT::RHSDriver */
	dArray2DT fIPVariable;
	/*@}*/

	//##########################################################################################
	//############## Attributes from ContinuumElementT.h needed for cut and paste ##############
	//############## methods in Traction_and_Body_Force.cpp (i.e. methods now in this class) ### 
	//##########################################################################################

	public:

		enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2  /**< diagonally lumped mass */ };

	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;

	protected:

	 	/** apply traction boundary conditions to the coarse scale equations */
		void ApplyTractionBC(void);

		/** add contribution from the body force */
		void AddBodyForce(LocalArrayT& body_force) const;
	
		/** element body force contribution 
	 * \param mass_type mass matrix type of ContinuumElementT::MassTypeT
	 * \param constM pre-factor for the element integral
	 * \param nodal nodal values. Pass NULL for no nodal values: [nen] x [ndof]
	 * \param ip_values integration point source terms. Pass NULL for no integration
	 *        point values : [nip] x [ndof] */
	void FormMa(MassTypeT mass_type, double constM, const LocalArrayT* nodal_values, const dArray2DT* ip_values);
	 		
	void EchoTractionBC(ifstreamT& in, ostream& out);
	// could also break up. Input and defaults(per output format) are
	// shared but the output of what each code means is class-dependent
	void EchoBodyForce(ifstreamT& in, ostream& out);

  /** update traction BC data for the coarse scale equations */
	void SetTractionBC(void);

	/* body force vector */
	const ScheduleT* fBodySchedule; /**< body force schedule */
	dArrayT fBody; /**< body force vector   */

	/* traction data */
	ArrayT<Traction_CardT> fTractionList;
	int fTractionBCSet;

	dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
	
};

inline const ShapeFunctionT& StaggeredMultiScaleT::ShapeFunction(void) const 
{
#if __option(extended_errorcheck)
	if (!fShapes)
		ExceptionT::GeneralFail("StaggeredMultiScaleT::ShapeFunction", "no shape functions");
#endif
	return *fShapes;
}


} // namespace Tahoe 
#endif /* _STAGGERED_MULTISCALE_T_H_ */



