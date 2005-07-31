/* $Id: ContinuumElementT.h,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
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

class ContinuumElementT: public ElementBaseT
{
public:

	/* constructor */
	ContinuumElementT(FEManagerT& fe_manager);

	/* destructor */
	virtual ~ContinuumElementT(void);
	
	/* accessors */
	int NumIP(void) const;
	const int& CurrIP(void) const;
	const ShapeFunctionT& ShapeFunction(void) const;	
	const LocalArrayT& InitialCoordinates() const;
	const LocalArrayT& Displacements() const;

	/* allocates space and reads connectivity data */
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
	virtual void SideSetToFacets(int block_ID, const iArray2DT& sideset,
		iArray2DT& facets);

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

	/* driver for nodal value calculations */
	virtual void ComputeNodalValues(const iArrayT& output) = 0;

private:

	/* construct output labels array */
	virtual void SetOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) = 0;
	virtual void GenerateOutputLabels(const iArrayT& counts,
		ArrayT<StringT>& labels) const = 0;

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

	/* material data */
	MaterialListT* fMaterialList; 	

	/* output control */
	iArrayT	fOutputCodes;
	  	
	/* body force vector */
	int	    fBodyForceLTf;
	dArrayT fBody;	  	

	/* traction data */
	ArrayT<Traction_CardT> fTractionList;
	int fTractionBCSet;

	/* shape functions */
	ShapeFunctionT* fShapes;
	
	/* arrays with local ordering */
	LocalArrayT	fLocInitCoords;	// initial coords with local ordering
	LocalArrayT fLocDisp;	    // displacements with local ordering
	
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
	if (!fShapes) throw eGeneralFail;
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
