/* $Id: APS_AssemblyT.h,v 1.6 2003-09-21 22:14:36 raregue Exp $ */ 
//DEVELOPMENT
#ifndef _APS_ASSEMBLY_T_H_ 
#define _APS_ASSEMBLY_T_H_ 

#include "ContinuumT.h"

/* base classes */
#include "ElementBaseT.h"
#include "StringT.h"
#include "ofstreamT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

/* base multiscale classes */
#include "APS_FEA.h"
#include "APS_EnumT.h"
#include "APS_VariableT.h"
#include "APS_Bal_EqT.h"
#include "APS_BCJT.h"
//#include "BalLinMomT.h"
//#include "PlastT.h"
#include "FEA_FormatT.h"

namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class Traction_CardT;	/** mass types */

/** APS_AssemblyT: This class contains methods pertaining to kinematics of
 * a dual field formulation for strict anti-plane shear. These include a scalar 
 * out-of-plane displacement u and plastic gradient gamma_p
 **/

class APS_AssemblyT: public ElementBaseT
{
 //----- Class Methods -----------
	
 public:

	enum fMat_T 	{ 
									k__mu,
									k__m_rate,
									k__l,
									k__H,
									k__gamma0_dot,
									k__m1,
									k__m2,
									kNUM_FMAT_TERMS	}; // MAT for material here, not matrix

	/** constructor */
	APS_AssemblyT(const ElementSupportT& support, const FieldT& displ, 
		const FieldT& gammap);

	/** destructor */
	~APS_AssemblyT(void);

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

	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;
	
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
private:

	/** Data at time steps n and n+1 used by both Coarse and Fine */
	//APS_VariableT n,np1; // <-- keep local scope in elmt loop for now 

	/** Gradients with respect to reference coodinates */
	FEA_dMatrixT fgrad_gamma_p, fgrad_gamma_p_n, fVars_matrix;
	FEA_dVectorT fgrad_u, fgrad_u_n, fgamma_p, fgamma_p_n, fVars_vector;

	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode;

	/** number of integration points */
	int	fNumIP;
	/*@}*/

	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT u;		//total out-of-plane displacement
	LocalArrayT u_n; 	//total out-of-plane displacement from previous increment
	LocalArrayT del_u;	//the Newton-R update i.e. del_u = u - u_n (u_{n+1}^{k+1} implied
	LocalArrayT DDu;    //coarse scale acceleration (used for body force)
	LocalArrayT gamma_p;		//plastic gradient
	LocalArrayT gamma_p_n;
	LocalArrayT del_gamma_p;	//the Newton-R update
	dArrayT		del_u_vec;  		// vector form 
	dArrayT		del_gamma_p_vec;	// vector form
	/*@}*/

	int n_ip, n_sd, n_df, n_en, n_en_x_n_df; 
	int n_np, n_el, n_comps;
	//int step_number_last_iter;
	//bool New_Step;
	int step_number;
	int iPlastModelType;
	
	dArrayT fForces_at_Node;
	bool bStep_Complete;
 	double time;
 	
 	void Get_Fd_ext 	( dArrayT &fFd_ext );
	

	//-- Material Parameters 

	dArrayT fMaterial_Data;
	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the 
	 * coordinates in APS_AssemblyT::fCurrCoords, which are the
	 * current coordinates */
	ShapeFunctionT* fShapes;
	
	FEA_ShapeFunctionT fFEA_Shapes;

	/** reference coordinates */
	LocalArrayT fInitCoords;     

	/** current coordinates */
	LocalArrayT fCurrCoords;
	/*@}*/

	/* Data Storage */
	ElementMatrixT fKdd, fKdeps;
	ElementMatrixT fKepsd, fKepseps;
	dArrayT 	fFd_int;
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
	const FieldT& fDispl;
	
	/** the plastic gradient field */
	const FieldT& fPlast;	

	/** equations per element for gradient plasticity. The balance equation is
	 * in ElementBaseT::fEqnos and are handled by ElementBaseT. */
	iArray2DT fEqnos_plast;

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
	 * APS_AssemblyT::RHSDriver */
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
	void FormMa(MassTypeT mass_type, double constM, const LocalArrayT* nodal_values, 
				const dArray2DT* ip_values);
	 		
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

inline const ShapeFunctionT& APS_AssemblyT::ShapeFunction(void) const 
{
#if __option(extended_errorcheck)
	if (!fShapes)
		ExceptionT::GeneralFail("APS_AssemblyT::ShapeFunction", "no shape functions");
#endif
	return *fShapes;
}


} // namespace Tahoe 
#endif /* _APS_ASSEMBLY_T_H_ */



