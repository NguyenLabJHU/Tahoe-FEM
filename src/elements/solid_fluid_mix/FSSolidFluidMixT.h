/* $Id: FSSolidFluidMixT.h,v 1.7 2006-10-25 23:34:19 ebrahimi Exp $ */ 
//DEVELOPMENT
#ifndef _FS_SOLID_FLUID_MIX_T_H_ 
#define _FS_SOLID_FLUID_MIX_T_H_ 

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

/** FSSolidFluidMixT: This class contains a coupled finite deformation solid fluid
 * Total Lagrangian formulation in 3D.  It is assumed the mixture is saturated with 
 * the fluid phase, and that the solid and fluid phases are incompressible, whereas
 * the mixture is not.  Currently, the implementation is limited to a simple 
 * hyper-viscoelastic constitutive model.
 **/

class FSSolidFluidMixT: public ElementBaseT
{
	
public:

	enum fMaterial_T 	{ 
	    kMu,
	    kLambda,
	    kAlpha,
	    kRho_sR0,
	    kRho_fR0,
	    kPhi_s0,
	    kPhi_f0,
	    kKf,
	    kNUM_FMATERIAL_TERMS	};
									
	enum fIntegrate_T 	{ 
	    kBeta,
	    kGamma,
	    kNUM_FINTEGRATE_TERMS	};								

	/** constructor */
	FSSolidFluidMixT(	const ElementSupportT& support );				

	/** destructor */
	~FSSolidFluidMixT(void);
	
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
				AutoArrayT<const RaggedArray2DT<int>*>& eq_theta);

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

	void Select_Equations ( const int &iBalLinMom, const int &iBalMass );

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
	dArrayT fgrad_theta, fgrad_theta_n;
	
	dMatrixT fShapeSolid, fShapeSolidGrad, fShapeSolidGrad_temp;
	dArrayT fShapeFluid;
	dMatrixT fShapeFluidGrad;
	
	dMatrixT fDefGrad, fDefGradInv, fDefGradInvMatrix;

	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode_displ, fGeometryCodeSurf_displ, 
	    fGeometryCode_press, fGeometryCodeSurf_press;
	int fGeometryCode_displ_int, fGeometryCodeSurf_displ_int;

	/** number of integration points */
	int	fNumIP_displ, fNumIPSurf_displ, fNumIP_press, fNumIPSurf_press, 
	    knum_d_state, knum_i_state, knumstress, knumstrain, num_sidesets;
	/*@}*/

	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT u;		//solid displacement
	LocalArrayT u_n; 	//solid displacement from time t_n
	LocalArrayT del_u;	//displacement increment including Newton-R update, i.e. del_u = u_{n+1}^{k+1} - u_n
	LocalArrayT press;	//fluid pore pressure
	LocalArrayT press_n;	//fluid pore pressure from time t_n
	LocalArrayT del_press;	//pore pressure increment
	dArrayT		del_u_vec;  	// vector form 
	dArrayT		del_press_vec;	// vector form
	/*@}*/

	// problem size definitions
	int n_en_displ, n_en_press, n_en_displ_x_n_sd, n_sd_x_n_sd;
	int n_el, n_sd, n_sd_surf, n_en_surf;
	
	int step_number;
	int iConstitutiveModelType;
	
	//name of output vector
	StringT output;
	
	dArrayT fForces_at_Node;
	bool bStep_Complete;
 	double time, kappa;
 	
 	void Get_Fd_ext 	( dArrayT &fFd_ext );
	
	//-- Material Parameters 
	dArrayT fMaterial_Params;
	
	//-- Newmark Time Integration Parameters 
	dArrayT fIntegration_Params;
	
	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the 
	 * coordinates in FSSolidFluidMixT::fCurrCoords, which are the
	 * current coordinates */
	ShapeFunctionT* fShapes_displ;
	ShapeFunctionT* fShapes_press;

	/** reference coordinates */
	LocalArrayT fInitCoords_displ, fInitCoords_press;     
	/** current coordinates */
	LocalArrayT fCurrCoords_displ, fCurrCoords_press;
	/*@}*/

	/* Data Storage */
	ElementMatrixT fKdd, fKdtheta;
	ElementMatrixT fKthetad, fKthetatheta;
	dArrayT 	fFd_int;
	dArrayT 	fFd_ext;
	dArrayT		fFtheta_int;
	dArrayT		fFtheta_ext;

	dArrayT		fGRAD_disp;
	dArrayT 	fDefGradInv_Vector;
	dArrayT 	fEffective_Kirchhoff_vector;
	
	dMatrixT	fDeformation_Gradient;
	dMatrixT	fCauchy_Green_tensor;
	dMatrixT	fCauchy_Green_tensor_Inverse;
	dMatrixT	fDeformation_Gradient_Inverse;
	dMatrixT	fDeformation_Gradient_Transpose;
	dMatrixT	fDefGradInv_grad_GRAD;
	dMatrixT	fDefGradInv_grad_GRAD_Transpose;
	dMatrixT	fIdentity_matrix;
	dMatrixT	fTest_matrix_A;
	dMatrixT	fTest_matrix_B;
	dMatrixT	fTest_matrix_C;
        dMatrixT	fEffective_Second_Piola_tensor;
        dMatrixT	fTemp_matrix;
        dMatrixT	fEffective_Kirchhoff_tensor;
        dMatrixT	fIota_temp_matrix;
        dMatrixT	fVarpi_temp_matrix;


	/** the solid displacement field */
	const FieldT* fDispl;
	
	/** the fluid pore pressure field */
	const FieldT* fPress;

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
	ArrayT<const iArray2DT*> fConnectivities_press;
	ArrayT<iArray2DT> fConnectivities_reduced;
	/*@}*/

	/** \name equations */
	/*@{*/
	ArrayT<iArray2DT> fEqnos_displ;
	ArrayT<iArray2DT> fEqnos_press;
	/*@}*/

	/** \name element cards */
	/*@{*/
	AutoArrayT<ElementCardT> fElementCards_displ;
	AutoArrayT<ElementCardT> fElementCards_press;
	/*@}*/

	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;

	/** integration point stresses. Calculated and stored during 
	 * FSSolidFluidMixT::RHSDriver */
	dArray2DT fIPVariable;
	/*@}*/

	/** \name prescribed plastic gradient side set ID */
	/*@{*/
	ArrayT<StringT> fSideSetID;
	
	/** prescribed pore pressure weight over the side set;
	    the direction is defined by {n1,n2,n3} ?? */
	ArrayT<double> fPorePressureWght;

	/** for each side set, the global nodes on the faces in the set */
	ArrayT<iArray2DT> fPorePressureFaces;
	
	/** equation numbers for the nodes on each face */ 
	ArrayT<iArray2DT> fPorePressureFaceEqnos;
	
	/** side set elements */ 
	ArrayT<iArrayT> fSideSetElements;

	/** side set faces */ 
	ArrayT<iArrayT> fSideSetFaces;
	/*@}*/
	
	/** write output for debugging */
	/*@{*/
	/** output file stream */
	ofstreamT fs_mix_out;
	
	/** line output formating variables */
	int outputPrecision, outputFileWidth;
	/*@}*/

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

inline const ShapeFunctionT& FSSolidFluidMixT::ShapeFunction(void) const 
{
#if __option(extended_errorcheck)
	if (!fShapes_displ)
	    ExceptionT::GeneralFail("FSSolidFluidMixT::ShapeFunction", "no displ shape functions");
	if (!fShapes_press)
	    ExceptionT::GeneralFail("FSSolidFluidMixT::ShapeFunction", "no press shape functions");
#endif
	return *fShapes_displ;
	return *fShapes_press;
}

} // namespace Tahoe 
#endif /* _FS_SOLID_FLUID_MIX_T_H_ */



