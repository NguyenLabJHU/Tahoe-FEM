/* $Id: ContinuumElementT.h,v 1.3.2.2 2001-06-28 01:24:11 paklein Exp $ */
/* created: paklein (10/22/1996)                                          */
/* Interface for a general continuum element type, meaning the presence   */
/* of shape functions, and the implied presence of a continuum mechanics  */
/* object, although no member exists since the continuum specifics are too */
/* dependent on derived class specifics.                                  */

#ifndef _CONTINUUM_ELEMENT_T_H_
#define _CONTINUUM_ELEMENT_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

/* forward declarations */
class MaterialListT;
class ShapeFunctionT;
class Traction_CardT;
class StringT;

/** base class for elements using shape functions */
class ContinuumElementT: public ElementBaseT
{
public:
	/** constructor.
	 * \param fe_manager used for system parameters */
	ContinuumElementT(FEManagerT& fe_manager);

	/** destructor */
	virtual ~ContinuumElementT(void);
		
	/** number of element integration points */
	int NumIP(void) const;
	
	/** reference to element shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;	

	/** reference to the current integration point number */
	const int& CurrIP(void) const;
	
	/** the coordinates of the current integration point */
	void IP_Coords(dArrayT& ip_coords) const;

	/** field gradients.
	 * compute the gradient of the field at the current integration point 
	 * \param field nodal values of the field 
	 * \param gradient field gradient: [ndof] x [nsd] */
	void IP_ComputeGradient(const LocalArrayT& field, dMatrixT& gradient) const;

	/** element coordinates.
	 * \return initial nodal coordinates of current element: [nen] x [nsd] */
	const LocalArrayT& InitialCoordinates() const;
	
	/** element displacements.
	 * \return nodal displacements of current element: [nen] x [ndof] */
	const LocalArrayT& Displacements() const;

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/* set element group for new global equations numbers */
	virtual void Reinitialize(void);

	/* form of tangent matrix - symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual void ResetStep(void); // restore last converged state

	/* restart operations */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/* side set to nodes on facets data - dimensions facets */
	void SideSetToFacets(int block_ID, const iArray2DT& sideset, iArray2DT& facets) const;

	/* return geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, iArrayT& num_facet_nodes) const;
	
	/* return the geometry code */
	GeometryT::CodeT GeometryCode(void) const;

	/* initial condition/restart functions (per time sequence) */
	virtual void InitialCondition(void);

	/* surface facets */
	void SurfaceFacets(GeometryT::CodeT& geometry,
		iArray2DT& surface_facets, iArrayT& surface_nodes) const;
	void SurfaceFacets(GeometryT::CodeT& geometry,
		ArrayT<iArray2DT>& surface_facet_sets,
		iArrayT& surface_nodes) const;
		//sorts facets into connected sets
	
	/* surface nodes */
	void SurfaceNodes(iArrayT& surface_nodes) const;
	
	/* materials list */
	const MaterialListT& MaterialsList(void) const;
	
protected:

	/* mass types */
	enum MassTypeT {kNoMass = 0,
            kConsistentMass = 1,
                kLumpedMass = 2};

	/* initialization functions */
	virtual void SetLocalArrays(void);
	virtual void SetShape(void) = 0;

	/* form the residual force vector */
	virtual void RHSDriver(void);

	/* compute contribution to RHS from traction BC's */
	void ApplyTractionBC(void);

	/*  form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/* accumulate the element mass matrix */
	void FormMass(int mass_type, double constM);

	/* add contribution from the body force */
	void AddBodyForce(LocalArrayT& body_force) const;
	
	/* element body force contribution */
	void FormMa(int mass_type, double constM, const LocalArrayT& body_force);
	 		
	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);	
	virtual void WriteMaterialData(ostream& out) const;
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out) = 0;
	void EchoBodyForce(ifstreamT& in, ostream& out);
	// could also break up. Input and defaults(per output format) are
	// shared but the output of what each code means is class-dependent
	void EchoTractionBC(ifstreamT& in, ostream& out);

	/* return a pointer to a new material list */
	virtual MaterialListT* NewMaterialList(int size) const = 0;

	/* return the "bounding" elements and the corresponding
	 * neighbors, both dimensioned internally */
	void BoundingElements(iArrayT& elements, iArray2DT& neighbors) const;

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

	/* driver for calculating output values */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values) = 0;

	/* check material outputs - return true if OK */
	virtual bool CheckMaterialOutput(void) const;

private:

	/* construct output labels array */
	virtual void GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const = 0;

	void EchoTractionBC_ASCII(ifstreamT& in, ostream& out);
	void EchoTractionBC_TahoeII(ifstreamT& in, ostream& out);
	void EchoTractionBC_ExodusII(ifstreamT& in, ostream& out);

	/* update traction BC data */
	void SetTractionBC(void);

	/* return the default number of element nodes */
	virtual int DefaultNumElemNodes(void) const;
	//NOTE: needed because ExodusII does not store ANY information about
	//      empty element groups, which causes trouble for parallel execution
	//      when a partition contains no element from a group.

protected:

	/* control data */
	GeometryT::CodeT fGeometryCode;
	int	fNumIP;
	int fOutputID;

	/* materials */
	MaterialListT*        fMaterialList;  /**< list of materials */
	ArrayT<ArrayT<bool> > fMaterialNeeds; /**< list of needs by material */
	
	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	  	
	/* body force vector */
	int	    fBodyForceLTf; /**< body force schedule */
	dArrayT fBody;	  	   /**< body force vector   */

	/* traction data */
	ArrayT<Traction_CardT> fTractionList;
	int fTractionBCSet;

	/** shape functions */
	ShapeFunctionT* fShapes;
	
	/* arrays with local ordering */
	LocalArrayT	fLocInitCoords;	/**< initial coords with local ordering */
	LocalArrayT fLocDisp;	    /**< displacements with local ordering  */ 
	
	/* work space */
	dArrayT fNEEvec; // [element DOF]
	dArrayT fDOFvec; // [nodal DOF]
	dArrayT fNSDvec; // [nodal dim]
};

/* inlines */

/* return the geometry code */
inline GeometryT::CodeT ContinuumElementT::GeometryCode(void) const
{ return fGeometryCode; }

/* accessors */
inline int ContinuumElementT::NumIP(void) const { return fNumIP; }

inline const ShapeFunctionT& ContinuumElementT::ShapeFunction(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes)
	{
		cout << "\n ContinuumElementT::ShapeFunction: no shape functions" << endl;
		throw eGeneralFail;
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
	if (!fMaterialList) throw eGeneralFail;
	return *fMaterialList;
}

#endif /* _CONTINUUM_ELEMENT_T_H_ */
