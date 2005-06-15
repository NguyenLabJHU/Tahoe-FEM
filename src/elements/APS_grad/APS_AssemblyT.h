/* $Id: APS_AssemblyT.h,v 1.33 2005-06-15 00:06:15 raregue Exp $ */ 
//DEVELOPMENT
#ifndef _APS_ASSEMBLY_T_H_ 
#define _APS_ASSEMBLY_T_H_ 

#include "ContinuumT.h"
#include "ModelManagerT.h"

/* base classes */
#include "ElementBaseT.h"
#include "StringT.h"
#include "Traction_CardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"

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

/* base multiscale classes */
#include "APS_FEA.h"
#include "APS_EnumT.h"
#include "APS_VariableT.h"
#include "APS_Bal_EqT.h"
#include "APS_kappa_alpha_macT.h"
#include "APS_kappa_alphaT.h"
#include "APS_kappa_gpT.h"
#include "FEA_FormatT.h"

namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class Traction_CardT;	
class StringT;

/** APS_AssemblyT: This class contains kinematics of
 * a dual field formulation for strict anti-plane shear. These include a scalar 
 * out-of-plane displacement u and an in-plane plastic gradient vector gamma_p.
 * This problem is two-dimensional.
 **/

class APS_AssemblyT: public ElementBaseT
{
//----- Class Methods -----------
	
public:

	enum fMat_T 	{ 
									kMu,
									km_rate,
									kl,
									kH,
									kgamma0_dot_1,
									kgamma0_dot_2,
									kgamma0_dot_3,
									km1_x,
									km1_y,
									km2_x,
									km2_y,
									km3_x,
									km3_y,									
									kkappa0_1,
									kkappa0_2,
									kkappa0_3,
									kNUM_FMAT_TERMS	}; // MAT for material here, not matrix

	/** constructor */
	APS_AssemblyT(	const ElementSupportT& support );				

	/** destructor */
	~APS_AssemblyT(void);
	
	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;

	/** echo input */
	void Echo_Input_Data (void); 

	/** return true if the element contributes to the solution of the
	 * given group. ElementBaseT::InGroup returns true if group is the
	 * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
	virtual bool InGroup(int group) const;

	/* initialize/finalize time increment */
	//virtual void InitStep(void);
	virtual void CloseStep(void);
	//virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state

	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations( AutoArrayT<const iArray2DT*>& eq_d,
					AutoArrayT<const RaggedArray2DT<int>*>& eq_eps);

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

	void Select_Equations ( const int &iBalLinMom, const int &iPlast );

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

	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	
private:

	/** Data at time steps n and n+1 used by both Coarse and Fine */
	//APS_VariableT n,np1; // <-- keep local scope in elmt loop for now 

	/** Gradients and other matrices */
	FEA_dMatrixT fgrad_gamma_p, fgrad_gamma_p_n, fVars_matrix, fgrad_u, fgrad_u_n, 
				fgrad_u_surf, fgrad_u_surf_n;
	FEA_dVectorT fgamma_p, fgamma_p_n, fgamma_p_surf, fgamma_p_surf_n, fVars_vector, fstate, fstate_n;

	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode_displ, fGeometryCodeSurf_displ, 
			fGeometryCode_plast, fGeometryCodeSurf_plast;
	int fGeometryCode_displ_int, fGeometryCodeSurf_displ_int;

	/** number of integration points */
	int	fNumIP_displ, fNumIPSurf_displ, fNumIP_plast, fNumIPSurf_plast, 
		knum_d_state, knum_i_state, knumstress, knumstrain, num_sidesets;
	/*@}*/

	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT u;		//total out-of-plane displacement
	LocalArrayT u_n; 	//total out-of-plane displacement from previous increment
	LocalArrayT del_u;	//the Newton-R update i.e. del_u = u - u_n (u_{n+1}^{k+1} implied)
	LocalArrayT gamma_p;		//plastic gradient
	LocalArrayT gamma_p_n;
	LocalArrayT del_gamma_p;	//the Newton-R update
	dArrayT		del_u_vec;  		// vector form 
	dArrayT		del_gamma_p_vec;	// vector form
	/*@}*/

	//int n_ip_displ, n_ip_plast;
	int n_en_displ, n_en_plast, n_en_plast_x_n_sd;
	int n_el, n_sd, n_sd_surf, n_en_surf;
	//int step_number_last_iter;
	//bool New_Step;
	int step_number;
	int iPlastModelType;
	
	//name of output vector
	StringT output;
	
	ArrayT < FEA_dVector_ArrayT >  Render_Vector; // ( n_el x num variables to render (n_rv)

	dArrayT fForces_at_Node;
	bool bStep_Complete;
 	double time, kappa;
 	
 	void Get_Fd_ext 	( dArrayT &fFd_ext );
	

	//-- Material Parameters 
	dArrayT fMaterial_Data;
	
	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the 
	 * coordinates in APS_AssemblyT::fCurrCoords, which are the
	 * current coordinates */
	ShapeFunctionT* fShapes_displ;
	ShapeFunctionT* fShapes_plast;
	
	// shape functions in FEA class type
	FEA_ShapeFunctionT fFEA_Shapes_displ, fFEA_Shapes_plast;
	FEA_SurfShapeFunctionT fFEA_SurfShapes;

	// reference and current coordinates are the same for anti-plane shear
	/** reference coordinates */
	LocalArrayT fInitCoords_displ, fInitCoords_plast;     
	/** current coordinates */
	LocalArrayT fCurrCoords_displ, fCurrCoords_plast;
	/*@}*/

	/* Data Storage */
	ElementMatrixT fKdd, fKdeps, fKdd_face;
	ElementMatrixT fKepsd, fKepseps;
	dArrayT 	fFd_int, fFd_int_face;
	dArrayT 	fFd_ext;
	dArrayT		fFeps_int;
	dArrayT		fFeps_ext;

	/* Multi-Field Element Formulators */
	BalLinMomT* fEquation_d;	
	PlastT* 	fEquation_eps;

	/* Multi-Field Materials */
	APS_MaterialT* fPlastMaterial;
	APS_MaterialT* fBalLinMomMaterial;

	/* Conversion methods: puts data in FEA format (very little cost in perspective) */
	FEA_FormatT Convert;

	/** the displacement field */
	const FieldT* fDispl;
	
	/** the plastic gradient field */
	const FieldT* fPlast;

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
	ArrayT<const iArray2DT*> fConnectivities_plast;
	ArrayT<iArray2DT> fConnectivities_reduced;
	/*@}*/

	/** \name equations */
	/*@{*/
	ArrayT<iArray2DT> fEqnos_displ;
	ArrayT<iArray2DT> fEqnos_plast;
	/*@}*/

	/** \name element cards */
	/*@{*/
	AutoArrayT<ElementCardT> fElementCards_displ;
	AutoArrayT<ElementCardT> fElementCards_plast;
	/*@}*/

	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;

	/** integration point stresses. Calculated and stored during 
	 * APS_AssemblyT::RHSDriver */
	dArray2DT fIPVariable;
	/*@}*/

	/** \name prescribed plastic gradient side set ID */
	/*@{*/
	ArrayT<StringT> fSideSetID;
	
	/** prescribed plastic gradient scalar weight over the side set;
	    the direction is defined by {m1,m2} */
	ArrayT<double> fPlasticGradientWght;

	/** for each side set, the global nodes on the faces in the set */
	ArrayT<iArray2DT> fPlasticGradientFaces;
	
	/** equation numbers for the nodes on each face */ 
	ArrayT<iArray2DT> fPlasticGradientFaceEqnos;
	
	/** side set elements */ 
	ArrayT<iArrayT> fSideSetElements;

	/** side set faces */ 
	ArrayT<iArrayT> fSideSetFaces;
	/*@}*/

	//##########################################################################################
	//############## Attributes from ContinuumElementT.h needed for cut and paste ##############
	//############## methods in APS_Traction.cpp (i.e. methods now in this class) ### 
	//##########################################################################################

public:

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
	LocalArrayT fLocDisp;	      /**< displacements with local ordering  */ 
	/*@}*/

	/** \name work space */
	/*@{*/
	dArrayT fNEEvec; /**< work space vector: [element DOF] */
	dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
	/*@}*/
	
};

inline const ShapeFunctionT& APS_AssemblyT::ShapeFunction(void) const 
{
#if __option(extended_errorcheck)
	if (!fShapes_displ)
		ExceptionT::GeneralFail("APS_AssemblyT::ShapeFunction", "no displ shape functions");
	if (!fShapes_plast)
		ExceptionT::GeneralFail("APS_AssemblyT::ShapeFunction", "no plast shape functions");
#endif
	return *fShapes_displ;
	return *fShapes_plast;
}



} // namespace Tahoe 
#endif /* _APS_ASSEMBLY_T_H_ */



