/* $Id: ContinuumElementT.h,v 1.26.24.2 2005-02-24 01:14:17 thao Exp $ */
/* created: paklein (10/22/1996) */
#ifndef _CONTINUUM_ELEMENT_T_H_
#define _CONTINUUM_ELEMENT_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class MaterialListT;
class MaterialSupportT;
class ShapeFunctionT;
class Traction_CardT;
class StringT;

/** base class for elements using shape functions */
class ContinuumElementT: public ElementBaseT
{
public:

	/** constructor */
	ContinuumElementT(const ElementSupportT& support, const FieldT& field);
	ContinuumElementT(const ElementSupportT& support);

	/** destructor */
	virtual ~ContinuumElementT(void);
		
	/** number of element integration points */
	int NumIP(void) const { return fNumIP;} ;
	
	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;

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

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** collecting element group equation numbers. This call from the FEManagerT
	 * is a signal to the element group that the equation system is up to date
	 * for the current time increment. See ElementBaseT::Equations for more
	 * information. */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** form of tangent matrix - symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);
	
	 /** write restart information to stream */
	virtual void WriteRestart(ostream& out) const;

	/** register self for output */
	virtual void RegisterOutput(void);

	/** send output */
	virtual void WriteOutput(void);

	/** return geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, iArrayT& num_facet_nodes) const;
	
	/** return the geometry code */
	GeometryT::CodeT GeometryCode(void) const;

	/*set active elements*/
	virtual void SetStatus(const ArrayT<StatusT>& status);

	/** initial condition/restart functions (per time sequence) */
	virtual void InitialCondition(void);
	
	/** reference to the materials list */
	const MaterialListT& MaterialsList(void) const;
	
	/** mass types */
	enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2  /**< diagonally lumped mass */ };

protected:

	/** stream extraction operator */
	friend istream& operator>>(istream& in, ContinuumElementT::MassTypeT& type);

	/** echo element connectivity data. Calls the inherited ElementBaseT::ElementBaseT
	 * and then constructs the communicator for the processes with non-zero numbers
	 * of elements in this group */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** allocate and initialize shape function objects */
	virtual void SetShape(void) = 0;

	/** form the residual force vector. computes contribution from natural
	 * boundary conditions */
	virtual void RHSDriver(void);

	/** compute contribution to element residual force due to natural boundary 
	 * conditions */
	void ApplyTractionBC(void);
	virtual bool Axisymmetric(void) const { return false; };

	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** accumulate the element mass matrix */
	void FormMass(int mass_type, double constM);

	/** add contribution from the body force */
	void AddBodyForce(LocalArrayT& body_force) const;
	
	/** element body force contribution 
	 * \param mass_type mass matrix type of ContinuumElementT::MassTypeT
	 * \param constM pre-factor for the element integral
	 * \param nodal nodal values. Pass NULL for no nodal values: [nen] x [ndof]
	 * \param ip_values integration point source terms. Pass NULL for no integration
	 *        point values : [nip] x [ndof] */
	void FormMa(MassTypeT mass_type, double constM, 
		const LocalArrayT* nodal_values,
		const dArray2DT* ip_values);
	 		
	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;
	
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);	
	virtual void WriteMaterialData(ostream& out) const;
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out) = 0;
	void EchoBodyForce(ifstreamT& in, ostream& out);
	// could also break up. Input and defaults(per output format) are
	// shared but the output of what each code means is class-dependent
	void EchoTractionBC(ifstreamT& in, ostream& out);

	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param nsd number of spatial dimensions
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(int nsd, int size) = 0;

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

	/** \name calculating output values */
	/*@{*/
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values) = 0;
	/*@}*/

	/** check consistency of material outputs.
	 * \return true if output variables of all materials for the group matches */
	virtual bool CheckMaterialOutput(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** construct output labels array */
	virtual void GenerateOutputLabels(
		const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const = 0;

	/** update traction BC data */
	void SetTractionBC(void);

	/** return the default number of element nodes.
	 * \note needed because ExodusII does not store \a any information about
	 * empty element groups, which causes trouble for parallel execution
	 * when a partition contains no elements from a group. */
	virtual int DefaultNumElemNodes(void) const;

protected:

	/** communicator over processes with elements in this group */
	CommunicatorT* fGroupCommunicator;

	/** list of materials */
	MaterialListT* fMaterialList;
	
	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	  	
	/* body force vector */
	const ScheduleT* fBodySchedule; /**< body force schedule */
	dArrayT fBody; /**< body force vector   */

	/* traction data */
	ArrayT<Traction_CardT> fTractionList;
	int fTractionBCSet;

	/** shape functions */
	ShapeFunctionT* fShapes;
	
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

private:

	/** number of integration points */
	int	fNumIP;

	/** output ID */
	int fOutputID;

	/** element parameter */
	GeometryT::CodeT fGeometryCode;
	
	/** cached results from ContinuumElementT::Axisymmetric */
	bool fAxisymmetric;
};

/* inlines */
/* communicator over the group */
inline const CommunicatorT& ContinuumElementT::GroupCommunicator(void) const
{
#if __option(extended_errorcheck)
	if (!fGroupCommunicator)
		ExceptionT::GeneralFail("ContinuumElementT::GroupCommunicator", "pointer not set");
#endif
	return *fGroupCommunicator;
}

/* return the geometry code */
inline GeometryT::CodeT ContinuumElementT::GeometryCode(void) const
{ return fGeometryCode; }

/* accessors */
inline const ShapeFunctionT& ContinuumElementT::ShapeFunction(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes)
	{
		cout << "\n ContinuumElementT::ShapeFunction: no shape functions" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif
	return *fShapes;
}

inline const LocalArrayT& ContinuumElementT::InitialCoordinates() const
{	
	return fLocInitCoords;
}

inline const LocalArrayT& ContinuumElementT::Displacements() const
{
	return fLocDisp;
}

inline const MaterialListT& ContinuumElementT::MaterialsList(void) const
{
	if (!fMaterialList) throw ExceptionT::kGeneralFail;
	return *fMaterialList;
}

} // namespace Tahoe 
#endif /* _CONTINUUM_ELEMENT_T_H_ */
