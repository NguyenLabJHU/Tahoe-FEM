#ifndef _BRIDGING_SCALE_T_H_
#define _BRIDGING_SCALE_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class MaterialListT;
class ShapeFunctionT;
class Traction_CardT;
class StringT;

/** base class for elements using shape functions */
class BridgingScaleT: public ElementBaseT
{
public:

	/** constructor */
	BridgingScaleT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	virtual ~BridgingScaleT(void);
		
	/** number of element integration points */
	int NumIP(void) const { return fNumIP;} ;
	
	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;

	/** interpolate the nodal field values to the current integration point */
    void IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u) const;

	/** interpolate the nodal field values to the specified integration point */
    void IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u, int ip) const;

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

	/* form of tangent matrix - symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual void ResetStep(void); // restore last converged state

	/** register self for output */
	virtual void RegisterOutput(void);

	/** send output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/** return the geometry code */
	GeometryT::CodeT GeometryCode(void) const;

	/** reference to the materials list */
	const MaterialListT& MaterialsList(void) const;
	
	/** mass types */
	enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2  /**< diagonally lumped mass */ };

protected:

	/** stream extraction operator */
	friend istream& operator>>(istream& in, BridgingScaleT::MassTypeT& type);

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** allocate and initialize shape function objects */
	virtual void SetShape(void) = 0;

	/** form the residual force vector. computes contribution from natural
	 * boundary conditions */
	virtual void RHSDriver(void);

	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;
	
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);	
	virtual void WriteMaterialData(ostream& out) const;
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out) = 0;

	/** construct a new material list and return a pointer */
	virtual MaterialListT* NewMaterialList(int size) const = 0;

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

	/* driver for calculating output values */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values) = 0;


private:

	/* construct output labels array */
	virtual void GenerateOutputLabels(
		const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const = 0;

	/** return the default number of element nodes.
	 * \note needed because ExodusII does not store \a any information about
	 * empty element groups, which causes trouble for parallel execution
	 * when a partition contains no elements from a group. */
	virtual int DefaultNumElemNodes(void) const;

	/* bridging scale-related computational functions */
	
	/* computes error caused by projecting solution onto FEM basis space */
	void ComputeError(void);

	/* computes "fine scale" displacement, ie MD - overlap due to FEM */
	void ComputeFineScaleU(void);

protected:

	/* materials */
	MaterialListT* fMaterialList;  /**< list of materials */
	
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
	
	/* arrays with local ordering */
	LocalArrayT	fLocInitCoords;	/**< initial coords with local ordering */
	LocalArrayT fLocDisp;	    /**< displacements with local ordering  */ 
	
	/* work space */
	dArrayT fNEEvec; /**< work space vector: [element DOF] */
	dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
//	dArrayT fNSDvec; /**< work space vector: [nodal dim]   */

private:

	/** number of integration points */
	int	fNumIP;

	/** output ID */
	int fOutputID;

	/* control data */
	GeometryT::CodeT fGeometryCode;
};

/* inlines */

/* return the geometry code */
inline GeometryT::CodeT BridgingScaleT::GeometryCode(void) const
{ return fGeometryCode; }

/* accessors */
inline const ShapeFunctionT& BridgingScaleT::ShapeFunction(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes)
	{
		cout << "\n BridgingScaleT::ShapeFunction: no shape functions" << endl;
		throw eGeneralFail;
	}
#endif
	return *fShapes;
}

inline const LocalArrayT& BridgingScaleT::InitialCoordinates() const
{	
	return fLocInitCoords;
}

inline const LocalArrayT& BridgingScaleT::Displacements() const
{
	return fLocDisp;
}

inline const MaterialListT& BridgingScaleT::MaterialsList(void) const
{
	if (!fMaterialList) throw eGeneralFail;
	return *fMaterialList;
}

} // namespace Tahoe 
#endif /* _BRIDGING_SCALE_T_H_ */
