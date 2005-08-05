/* $Id: MFGPAssemblyT.h,v 1.2 2005-08-05 07:20:19 kyonten Exp $ */ 
//DEVELOPMENT
#ifndef _MFGP_ASSEMBLY_T_H_ 
#define _MFGP_ASSEMBLY_T_H_ 

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe 
{
/* forward declarations */
class Traction_CardT;
class StringT;
class D3MeshFreeShapeFunctionT;
class MFGPMatListT;
class MFGPMatSupportT;

/** MFGPAssemblyT: This class contains kinematics of
 * a dual field formulation for meshfree implementation of a gradient
 * plasticity model. These include a vector displacement u and scalar
 * plastic multiplier, lambda.
 * This model doesn't support 1D case.
 **/

//class MFGPAssemblyT: public MeshFreeFractureSupportT, public ElementBaseT
class MFGPAssemblyT: public ElementBaseT
{
	
public:
	/** constructor */
	MFGPAssemblyT(const ElementSupportT& support);

	/** destructor */
	virtual ~MFGPAssemblyT(void);

	/** return true if the element contributes to the solution of the
	 * given group. ElementBaseT::InGroup returns true if group is the
	 * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
	virtual bool InGroup(int group) const;
	
	/** number of element integration points */
	int NumIP(void) const { return fNumIP_displ;} ;
	
	/** reference to element shape functions */
	const D3MeshFreeShapeFunctionT& ShapeFunction(void) const;
	
	/** reference to the current integration point number */
	const int& CurrIP(void) const;

	/** communicator over the group */
	const CommunicatorT& GroupCommunicator(void) const;
	
	/** the coordinates of the current integration point */
	void IP_Coords(dArrayT& ip_coords) const;

	/** interpolate the nodal field values to the current integration point */
    void IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u) const;

	/** interpolate the nodal field values to the specified integration point */
    void IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u, int ip) const;

	/** field gradients.
	 * compute the gradient of the field at the current integration point 
	 * \param field nodal values of the field 
	 * \param gradient field gradient: [ndof] x [nsd] */
	void IP_ComputeGradient(const LocalArrayT& field, dMatrixT& gradient) const;

	/** field gradients.
	 * compute the gradient of the field at the specified integration point 
	 * \param field nodal values of the field 
	 * \param gradient field gradient: [ndof] x [nsd] */
	void IP_ComputeGradient(const LocalArrayT& field, dMatrixT& gradient, int ip) const;
	
	/** field Laplacian.
	 * compute the Laplacian of the field at the current integration point 
	 * \param field nodal values of the field 
	 * \param field Laplacian: [ndof] */
	void IP_ComputeLaplacian(const LocalArrayT& field, dArrayT& laplacian) const;
	
	/** field Laplacian.
	 * compute the Laplacian of the field at the specific integration point 
	 * \param field nodal values of the field 
	 * \param field Laplacian: [ndof] */
	void IP_ComputeLaplacian(const LocalArrayT& field, dArrayT& laplacian, int ip) const;
										
	/** extrapolate all integration point values to the nodes
	 * \param IPvalues values from the integration points: [nip] 
	 * \param nodalvalues extrapolated values: [nnd] */
	void IP_ExtrapolateAll(const dArrayT& ip_values, dArrayT& nodal_values) const;

	/** element coordinates.
	 * \return initial nodal coordinates of current element: [nen] x [nsd] */
	const LocalArrayT& InitialCoordinates() const;
	
	/** element displacements.
	 * \return nodal displacements of current element: [nen] x [ndof] */
	const LocalArrayT& Displacements() const;
	const LocalArrayT& LastDisplacements(void) const;
	const LocalArrayT& Accelerations(void) const;
	const LocalArrayT& PlasticMultipliers() const;
	const LocalArrayT& LastPlasticMultipliers(void) const;
	
	/** reference to the materials list */
	const MFGPMatListT& MFGPMatList(void) const;
	
	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations( AutoArrayT<const iArray2DT*>& eq_u,
							AutoArrayT<const RaggedArray2DT<int>*>& eq_lambda);
	
	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state
	
	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;
	
	/** register element for output */
	virtual void RegisterOutput(void);

	/** write element output */
	virtual void WriteOutput(void);	

	/** element level reconfiguration for the current time increment */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);
	/*@}*/
	
	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	
	/** return geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, iArrayT& num_facet_nodes) const;
	
	/** return the geometry code */
	GeometryT::CodeT GeometryCode(void) const;

	 /*set active elements*/
    virtual void SetStatus(const ArrayT<ElementCardT::StatusT>& status);
    
    /** initial condition/restart functions (per time sequence) */
	virtual void InitialCondition(void);
	
	/** mass types */
	/** mass types */
	enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2, /**< diagonally lumped mass */
             kAutomaticMass = 3  /**< select the mass type base on the time integration scheme */};
	MassTypeT static int2MassTypeT(int i);
	
	/* extrapolate the integration point stresses and internal variables, 
   	 * check the yield condition on the nodes of each background grid (element), 
     * and pass the flag whether the nodes are elastically or plastically loaded
	*/
	/*@{*/  
	virtual void CheckNodalYield(void);
	/*@}*/
	
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

	/** stream extraction operator */
	//friend istream& operator>>(istream& in, MFGPAssemblyT::MassTypeT& type);

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** allocate and initialize shape function objects */
	virtual void SetShape(void) = 0;
	
	/** compute contribution to element residual force due to natural boundary 
	 * conditions */
	void ApplyTractionBC(void);
	virtual bool Axisymmetric(void) const { return false; };
	
	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void);
	
	/** \name element loop operations */
	/*@{*/
	/** reset loop */
	virtual void Top(void);
	
	/* increment current element */
	virtual bool NextElement(void);
	
	/** form the element stiffness matrix
	 * Compute the linearization of the force calculated by MFGPAssemblyT::FormKd */
	virtual void FormStiffness(double constK) = 0;
	
	/** internal force */
	virtual void FormKd(double constK) = 0;

	/** accumulate the element mass matrix */
	void FormMass(int mass_type, double constM, bool axisymmetric);
	
	/** add contribution from the body force */
	void AddBodyForce(LocalArrayT& body_force) const;
	
	/** element body force contribution 
	 * \param mass_type mass matrix type of MFGPAssemblyT::MassTypeT
	 * \param constM pre-factor for the element integral
	 * \param nodal nodal values. Pass NULL for no nodal values: [nen] x [ndof]
	 * \param ip_values integration point source terms. Pass NULL for no integration
	 *        point values : [nip] x [ndof] */
	void FormMa(MassTypeT mass_type, double constM, bool axisymmetric,
		const LocalArrayT* nodal_values,
		const dArray2DT* ip_values);

	/* element data */
	void EchoBodyForce(ifstreamT& in, ostream& out);

	/** extract natural boundary condition information */
	void TakeNaturalBC(const ParameterListT& list);

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

	/** \name materials list methods
	 * These methods should be overridden to enable features derived from support
	 * for the MFGPMatListT data member MFGPAssemblyT::fMFGPMatList. To
	 * construct the material list, derived classes must override MFGPAssemblyT::CollectMaterialInfo
	 * to extract material parameters from the class parameters. These are then passed
	 * to MFGPAssemblyT::NewMFGPMatList to construct the materials list. */
	/*@{*/
	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;
	
	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param name list identifier
	 * \param size length of the list */
	virtual MFGPMatListT* NewMFGPMatList(const StringT& name, int size);
	
	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MFGPMatSupportT to be initialized. If NULL, allocate
	 *        a new MFGPMatSupportT and initialize it. */
	virtual MFGPMatSupportT* NewMFGPMatSupport(MFGPMatSupportT* p = NULL) const;
	
	/** check consistency of material outputs.
	 * \return true if output variables of all materials for the group matches */
	virtual bool CheckMaterialOutput(void) const;
	
	/** \name output methods 
	 * These methods should be overridden to enable output using the MFGPAssemblyT
	 * conventions. These methods are used during MFGPAssemblyT::RegisterOutput
	 * and MFGPAssemblyT::WriteOutput to send data for output. Two methods must
	 * be overridden to avoid execution of the methods inherited from ElementBaseT.
	 * -# MFGPAssemblyT::GenerateOutputLabels is used to create list of
	 *    output labels if given a non-empty list of output codes by MFGPAssemblyT::SetNodalOutputCodes
	 *    or MFGPAssemblyT::SetElementOutputCodes.
	 * -# MFGPAssemblyT::ComputeOutput is used to compute nodal and
	 *    element output data if given a non-empty list of output codes by MFGPAssemblyT::SetNodalOutputCodes
	 *    or MFGPAssemblyT::SetElementOutputCodes. */
	/*@{*/
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;

	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

	/** construct output labels array */
	virtual void GenerateOutputLabels(
		const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const;
	/*@}*/

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/
	
private:

	/** update traction BC data */
	void SetTractionBC(void);

	/** return the default number of element nodes.
	 * \note needed because ExodusII does not store \a any information about
	 * empty element groups, which causes trouble for parallel execution
	 * when a partition contains no elements from a group. */
	virtual int DefaultNumElemNodes(void) const;
	
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

	/** communicator over processes with elements in this group */
	CommunicatorT* fGroupCommunicator;
	
	/** list of materials */
	MFGPMatListT* fMFGPMatList;
	
	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	
	/* body force vector */
	const ScheduleT* fBodySchedule; /**< body force schedule */
	dArrayT fBody; /**< body force vector   */

	/* traction data */
	ArrayT<Traction_CardT> fTractionList;
	int fTractionBCSet;
	
	/** \name work space */
	/*@{*/
	dArrayT fNEEvec; /**< work space vector: [element DOF] */
	dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
	/*@}*/
	
	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;
	
	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT u;		//total displacement
	LocalArrayT u_n; 	//total displacement from previous increment
	LocalArrayT del_u;	//the Newton-R update i.e. del_u = u - u_n (u_{n+1}^{k+1} implied)
	LocalArrayT DDu;    //accelleration (used for body force)
	LocalArrayT lambda;		//plastic multiplier
	LocalArrayT lambda_n;
	LocalArrayT del_lambda;	//the Newton-R update
	dArrayT		del_u_vec;  	// vector form 
	dArrayT		del_lambda_vec;	// vector form
	/*@}*/
	
	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the 
	 * coordinates in MFGPAssemblyT::fCurrCoords, which are the
	 * current coordinates */
	 //different class for meshfree interpolations
	D3MeshFreeShapeFunctionT* fShapes_displ;
	D3MeshFreeShapeFunctionT* fShapes_plast;
	
	/** store shape function derivatives */
	bool fStoreShape;
	
	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode_displ, fGeometryCode_plast;
	
	// reference and current coordinates
	/** reference coordinates */
	LocalArrayT fInitCoords_displ, fInitCoords_plast; /* geometry */       
	/** current coordinates */
	LocalArrayT fCurrCoords_displ, fCurrCoords_plast;
	/*@}*/

	/** number of integration points */
	int	fNumIP_displ, fNumIP_plast; 
	
	double penalty_num; // to enforce b.c. on plastic multiplier (lambda)

	/* Data Storage */
	ElementMatrixT fKuu, fKuu_temp; // [ndof_displ]x[ndof_displ]; ndof_displ = nen_displ x dof_displ
	ElementMatrixT fKulambda, fKulambda_temp; // [ndof_displ]x[nen_plast]
	ElementMatrixT fKlambdau, fKlambdau_temp; // [nen_plast]x[ndof_displ]
	ElementMatrixT fKlambdalambda, fKlambdalambda_temp; // [nen_plast]x[nen_plast]; dof_plast = 1
	dArrayT 	fFu_int, fFu_int_temp; //[ndof_displ]
	dArrayT 	fFu_ext; //[ndof_displ]
	dArrayT		fFlambda, fFlambda_temp; // [nen_plast]	

	/** the displacement field */
	const FieldT* fDispl;
	
	/** the plastic mulitplier field */
	const FieldT* fPlast;	

	/** \name connectivities */
	/*@{*/
	ArrayT<const iArray2DT*> fConnectivities_displ;
	ArrayT<const iArray2DT*> fConnectivities_plast;
	/*@}*/

	/** \name element cards */
	/*@{*/
	AutoArrayT<ElementCardT> fElementCards_displ;
	AutoArrayT<ElementCardT> fElementCards_plast;
	/*@}*/

	/** integration point stresses. Calculated and stored during 
	 * MFGPAssemblyT::RHSDriver */
	dArray2DT fIPVariable;
	/*@}*/
	
	/* flags for nodal yield condition
	 * 1 if nodes satisfy yield condition */
	iArrayT fNodalYieldFlags;
	
	/* flags for penalty number 
	 * 1 if the penalty number is added */
	iArrayT fPenaltyFlags;
};

inline const D3MeshFreeShapeFunctionT& MFGPAssemblyT::ShapeFunction(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes_displ)
		ExceptionT::GeneralFail("MFGPAssemblyT::ShapeFunction", "no displ shape functions");
	if (!fShapes_plast)
		ExceptionT::GeneralFail("MFGPAssemblyT::ShapeFunction", "no plast shape functions");
#endif
		return *fShapes_displ;
		return *fShapes_plast;
}

/* inlines */
/* communicator over the group */
inline const CommunicatorT& MFGPAssemblyT::GroupCommunicator(void) const
{
#if __option(extended_errorcheck)
	if (!fGroupCommunicator)
		ExceptionT::GeneralFail("MFGPAssemblyT::GroupCommunicator", "pointer not set");
#endif
	return *fGroupCommunicator;
}

/* return the geometry code */
inline GeometryT::CodeT MFGPAssemblyT::GeometryCode(void) const
{return fGeometryCode_displ;}

inline const MFGPMatListT& MFGPAssemblyT::MFGPMatList(void) const
{
#if __option(extended_errorcheck)
	if (!fMFGPMatList) 
		ExceptionT::GeneralFail("MFGPAssemblyT::MFGPMatList", "no material list");
#endif
	return *fMFGPMatList;
}

inline const LocalArrayT& MFGPAssemblyT::InitialCoordinates() const {return fInitCoords_displ;}
inline const LocalArrayT& MFGPAssemblyT::Displacements() const {return u;}
inline const LocalArrayT& MFGPAssemblyT::LastDisplacements(void) const { return u_n; }
inline const LocalArrayT& MFGPAssemblyT::Accelerations(void) const { return DDu; }
inline const LocalArrayT& MFGPAssemblyT::PlasticMultipliers() const {return lambda;}
inline const LocalArrayT& MFGPAssemblyT::LastPlasticMultipliers(void) const { return lambda_n; }

} // namespace Tahoe 
#endif /* _MFGP_ASSEMBLY_T_H_ */

