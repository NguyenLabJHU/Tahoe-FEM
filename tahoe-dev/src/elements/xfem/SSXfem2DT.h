/* $Id: SSXfem2DT.h,v 1.1 2009-05-11 21:33:34 regueiro Exp $ */ 
//DEVELOPMENT
#ifndef _SS_XFEM_2D_T_H_ 
#define _SS_XFEM_2D_T_H_ 

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

/** SSXfem2DT: This class contains a coupled finite element 
 * implementation (displacement and enhanced displacement dofs)
 * in 2D for crack nucleation and propagation.  
 **/

class SSXfem2DT: public ElementBaseT
{
	
public:
/* material parameters */
	enum fMaterial_T 	{ 
		kMu,
		kLambda,
	    //add to this list
	    kNUM_FMATERIAL_TERMS	};
									
//	enum fIntegrate_T 	{ 
//	    kBeta,
//	    kGamma,
//	    kNUM_FINTEGRATE_TERMS	};								

	/** constructor */
	SSXfem2DT(	const ElementSupportT& support );				

	/** destructor */
	~SSXfem2DT(void);
	
	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunctionDispl(void) const;
	const ShapeFunctionT& ShapeFunctionEnhan(void) const;

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
				AutoArrayT<const RaggedArray2DT<int>*>& eq_q);

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

	void Select_Equations ( const int &iBalLinMom, const int &iBalLinMomEnh );

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
	dMatrixT fgrad_ujump, fgrad_ujump_n;
	
	dMatrixT fShapeDispl, fShapeDisplGrad;
	dMatrixT fShapeEnhan, fShapeEnhanGrad;
	
	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode_displ, fGeometryCodeSurf_displ, 
	    fGeometryCode_enhan, fGeometryCodeSurf_enhan;
	int fGeometryCode_displ_int, fGeometryCodeSurf_displ_int;

	/** number of integration points */
	int	fNumIP_displ, fNumIPSurf_displ, fNumIP_enhan, fNumIPSurf_enhan;
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
	LocalArrayT ujump;	//enhan-displacement-gradient dof
	LocalArrayT ujump_dot;	//enhan-displacement-gradient first time derivative
	LocalArrayT ujump_dot_n;	//enhan-displacement-gradient first time derivative from time t_n
	LocalArrayT ujump_dotdot;	//enhan-displacement-gradient second derivative  
	LocalArrayT ujump_dotdot_n;	//enhan-displacement-gradient second derivative from time t_n 
	LocalArrayT ujump_n;	//enhan-displacement-gradient from time t_n
	LocalArrayT del_ujump;	//enhan-displacement-gradient increment
	dArrayT		del_u_vec;  	// vector form 
	dArrayT		del_ujump_vec;	// vector form
	dArrayT		u_vec;  	// solid displacement in vector form 
	dArrayT		u_dot_vec;  	// solid velocity in vector form 
	dArrayT		u_dotdot_vec;  	// solid acceleration in vector form 
	dArrayT		ujump_vec;	// enhan-displacement-gradient in vector form
	dArrayT		ujump_dot_vec;	// first derivative of enhan-displacement-gradient in vector form
	dArrayT		ujump_dotdot_vec;	// second derivative of enhan-displacement-gradient in vector form

	/*@}*/

	// problem size definitions
	int n_en_displ, n_en_displ_x_n_sd, n_sd_x_n_sd;
	int n_el, n_sd, n_sd_surf, n_en_surf;
	int n_en_enhan, ndof_per_nd_enhan, n_en_enhan_x_ndof_per_nd_enhan, ndof_per_nd_enhan_x_n_sd;
	
	int step_number;
	int iConstitutiveModelType;
	
	//name of output vector
	StringT output;
	
	dArrayT fForces_at_Node;
	bool bStep_Complete;
	
 	double time;
 	
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
	ShapeFunctionT* fShapes_enhan;

	/** reference coordinates */
	LocalArrayT fInitCoords_displ, fInitCoords_enhan;     
	/** current coordinates */
	LocalArrayT fCurrCoords_displ, fCurrCoords_enhan;
	/*@}*/

	/* Data Storage */
	ElementMatrixT fKdd, fKdq;
	ElementMatrixT fKqd, fKqq;
	dArrayT 	fFd_int;
	dArrayT 	fFd_ext;
	dArrayT		fFq_int;
	dArrayT		fFq_ext;
	dArrayT		fTemp_vector_ndof_se;

	dMatrixT	fB_matrix;
	dMatrixT	fD_matrix;
	dMatrixT	fK_dd_BTDB_matrix;
	dMatrixT	fTemp_matrix_ndof_se_x_ndof_se;
	
	/* to store fstrain_vector_current_IP */
	dArrayT		fStrain_vector_current_IP;
	/* to store fstrain_IPs for each of the IPs of each element */
	dArray2DT	fStrain_IPs;
	/* to store fstress_vector_current_IP */
	dArrayT		fStress_vector_current_IP;
	/* to store fstress_IPs for each of the IPs of each element */
	dArray2DT	fStress_IPs;
	/* to store fState_IPs for each of the IPs of each element */
	dArray2DT	fState_variables_IPs;
	dArray2DT	fStrain_Elements_IPs;
	dArray2DT	fStress_Elements_IPs;
	dArray2DT	fState_variables_Elements_IPs;

	/** the solid displacement field */
	const FieldT* fDispl;
	
	/** the enhan-displacement-gradient field */
	const FieldT* fEnhan;

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
	ArrayT<const iArray2DT*> fConnectivities_enhan;
	ArrayT<iArray2DT> fConnectivities_reduced;
	/*@}*/

	/** \name equations */
	/*@{*/
	ArrayT<iArray2DT> fEqnos_displ;
	ArrayT<iArray2DT> fEqnos_enhan;
	/*@}*/

	/** \name element cards */
	/*@{*/
	AutoArrayT<ElementCardT> fElementCards_displ;
	AutoArrayT<ElementCardT> fElementCards_enhan;
	/*@}*/

	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;

	/** integration point stresses. Calculated and stored during 
	 * SSXfem2DT::RHSDriver */
	dArray2DT fIPVariable;
	/*@}*/

	/** \name prescribed plastic gradient side set ID */
	/*@{*/
	ArrayT<StringT> fSideSetID;
	
	/** prescribed enhan-displacement-gradient weight over the side set;
	    the direction is defined by {n1,n2,n3} ?? */
	//ArrayT<iArray2DT> fEnhanWght;

	/** for each side set, the global nodes on the faces in the set */
	ArrayT<iArray2DT> fEnhanFaces;
	
	/** equation numbers for the nodes on each face */ 
	ArrayT<iArray2DT> fEnhanFaceEqnos;
	
	/** side set elements */ 
	ArrayT<iArrayT> fSideSetElements;

	/** side set faces */ 
	ArrayT<iArrayT> fSideSetFaces;
	/*@}*/
	
	/** write output for debugging */
	/*@{*/
	/** output file stream */
	ofstreamT ss_xfem2D_out;
	
	/** line output formating variables */
	int outputPrecision, outputFileWidth;
	/*@}*/

	void Form_solid_shape_functions(const double* &shapes_displ_X);
	void Form_enhan_shape_functions(const double* &shapes_enhan_X);
	
	void Form_D_matrix(void);
	void Form_B_matrix(void);

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


inline const ShapeFunctionT& SSXfem2DT::ShapeFunctionDispl(void) const 
{
#if __option(extended_errorcheck)
	if (!fShapes_displ)
	    ExceptionT::GeneralFail("SSXfem2DT::ShapeFunctionDispl", "no displacement shape functions");
#endif
	return *fShapes_displ;
}

inline const ShapeFunctionT& SSXfem2DT::ShapeFunctionEnhan(void) const 
{
#if __option(extended_errorcheck)
	if (!fShapes_enhan)
	    ExceptionT::GeneralFail("SSXfem2DT::ShapeFunctionEnhan", "no enhanced shape functions");
#endif
	return *fShapes_enhan;
}

/* return the geometry code */
inline GeometryT::CodeT SSXfem2DT::GeometryCode(void) const
{ return fGeometryCode_displ; }


} // namespace Tahoe 
#endif /* _SS_XFEM_2D_T_H_ */



