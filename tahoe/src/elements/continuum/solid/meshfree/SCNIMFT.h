/* $Id: SCNIMFT.h,v 1.19 2004-10-26 22:07:51 paklein Exp $ */
#ifndef _SCNIMF_T_H_
#define _SCNIMF_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "VariArrayT.h"
#include "LocalArrayT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "MaterialListT.h"
#include "nVariArray2DT.h"
#include "InverseMapT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "GeometryBaseT.h"
#include "ScheduleT.h"
#include "LinkedListT.h"

#ifdef __QHULL__
#include "CompGeomT.h"
#endif

namespace Tahoe {

/** forward declarations */
class iGridManagerT;
class CommManagerT;
class dSPMatrixT; //TEMP
class InverseMapT;
class ifstreamT;
class ofstreamT;
class MeshFreeSupportT;
class Traction_CardT;

/** base class for particle types */
class SCNIMFT: public ElementBaseT
{
public:

	/** constructor */
	SCNIMFT(const ElementSupportT& support, const FieldT& field);
	SCNIMFT(const ElementSupportT& support);

	/** destructor */
	~SCNIMFT(void);
	
	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
			
	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) { return 0.0; };
	
	/** resgiter for writing output. Uses output labels generated by 
	 * small- and finite- strain implementations of this base class.
	 */
	virtual void RegisterOutput(void);

	/** write output. ParticleT::WriteOutput only writes search grid statistics.
	 * Sub-classes are responsible for writing data for each particle, given the
	 * variables names returned by ParticleT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

	/** Loop over nodes and compute stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type) = 0;
	
	/** Loop over nodes and compute internal force */
	virtual void RHSDriver(void);
	
	/** Generate local equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
						AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** Translate global node numbers to local ones -- communication routine for MFLagMultT */
	/** returns 0 if unsucessful, i.e. nodes not contained in fNodes */
	int GlobalToLocalNumbering(iArrayT& nodes);
	
	/* Translate global node numbers to local ones -- communication routine for MFLagMultT */
	int GlobalToLocalNumbering(RaggedArray2DT<int>& nodes);

	/** Return interpolated displacement field at selected nodes -- communication routine for MFLagMultT */
	void InterpolatedFieldAtNodes(const iArrayT& nodes, dArray2DT& fieldAtNodes);

	/** Return the data structure holding the supports of the localNodes and their window function values 
		-- communication routine for for MFLagMultT */
	void NodalSupportAndPhi(iArrayT& localNodes, RaggedArray2DT<int>& support, RaggedArray2DT<double>& phi);
	
	int SupportSize(int localNode);

	/** \name types needed for the Voronoi diagram calculation */
	/*@{*/
#ifndef __QHULL__
	/** Basic structure -- hullMap[i][j] gives index of jth vertex of structure i */
        typedef ArrayT<iArrayT> ConvexHullMap;

        /** Voronoi diagram facet structure is an array of convex hulls  */
        typedef ArrayT< ConvexHullMap > VoronoiDiagramMap;
#endif
	/*@}*/
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	virtual void TakeParameterList(const ParameterListT& list);
	
	void TakeNaturalBC(const ParameterListT& list);
	/*@}*/


protected: /* for derived classes only */

	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);
	
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const = 0;
	virtual MaterialListT* NewMaterialList(const StringT& name, int size) = 0;
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const;

	/** assemble particle mass matrix into LHS of global equation system */
	void AssembleParticleMass(const double rho);
	
	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(void);
	
	/** write out Voronoi diagram data */
	void VoronoiDiagramToFile(ofstreamT& vout);
	
	/** read in Voronoi diagram data */
	void VoronoiDiagramFromFile(ifstreamT& vin);
	
protected:

	MeshFreeSupportT* fMFSupport;
	
	/** pointer to list parameters needed to construct shape functions */
	const ParameterListT* fMeshfreeParameters;

	/** reference ID for sending output */
	int fOutputID;
	
	/** connectivities used to define the output set. */
	iArray2DT fPointConnectivities;

	/** \name cached RHS workspace */
	/*@{*/
	dArray2DT fForce;
	nVariArray2DT<double> fForce_man;
	/*@}*/
	
	/** spatial dimensionality */
	int fSD;

	/** indices of nodes */
	iArrayT fNodes;
	/** the coordinates of the nodes */
	dArray2DT fDeloneVertices;
	
	/** \name Geometrical Data Structures */
	/*@{*/
	
	/** these are dual to Voronoi facets. They have minor dimension of 2 . Difference in the two points is parallel to the normal vector of the dual facet. */
	iArray2DT fDeloneEdges;

	/** Voronoi facets dual to the Delone Edges */
	iArray2DT fDualFacets; // Tag for Deletion
	iArray2DT fSelfDualFacets; 

	/** Self-dual facet information. I.E. facets that contribute only to one integral over one boundary node's cell */
#ifdef __QHULL__
	CompGeomT::ConvexHullMap fSelfDuals;
#else
	ConvexHullMap fSelfDuals; // Tag for Deletion
#endif
	int fNumSelfDuals;
	int fNumClippedFacets;

	/** connectivity of boundary nodes. Currently determined from an underlying 
	    element connectivity */
	iArray2DT fBoundaryConnectivity; // Tag for Deletion
	
	/** union of nodes in fBoundaryConnectivity */
	iArrayT fBoundaryNodes; // Tag for Deletion
	
	/** true if boundary connectivity is simplicial */
	bool fBoundaryIsTriangulated; // Tag for Deletion
	
	/** additional edges associated only with one node */
	iArrayT fNonDeloneEdges; 
	
	/** normal vectors of the facets for those edges */
	dArray2DT fNonDeloneNormals;
	
	/** areas of boundary facets */
	dArrayT fBoundaryIntegrationWeights;

	/** Compute or read the Voronoi Diagram */	
#ifdef __QHULL__	
	CompGeomT* fVoronoi;
	CompGeomT::ConvexHullMap fVoronoiCells;
	CompGeomT::VoronoiDiagramMap fVoronoiFacetIndices;
#else
	void* fVoronoi;
	ConvexHullMap fVoronoiCells; // Tag for Deletion
	VoronoiDiagramMap fVoronoiFacetIndices; // Tag for Deletion
#endif
	bool qComputeVoronoiCell;
	StringT vCellFile;
	
	int fNumIP;

	ArrayT<dArrayT> fVoronoiFacetAreas; // Tag for Deletion
	ArrayT<dArray2DT> fVoronoiFacetNormals; // Tag for Deletion

	/** Volume associated with each node -- integration weight for nodal integration */
	dArrayT fVoronoiCellVolumes;
	dArray2DT fVoronoiVertices; // Tag for Deletion
	/*@}*/

	/** list of materials */
	MaterialListT* fMaterialList;
	
	/** workspaces for strain smoothing */
	ArrayT< LinkedListT<int> > nodeWorkSpace; // should be local?
	ArrayT< LinkedListT<dArrayT> > facetWorkSpace; // should be local?

	RaggedArray2DT<int> nodalCellSupports;
	RaggedArray2DT<dArrayT> bVectorArray;
	
	/** workspace for nodal shape functions */
	RaggedArray2DT<double> fNodalPhi, fBoundaryPhi;
	RaggedArray2DT<int> fNodalSupports, fBoundarySupports;
	  	
	/* body force vector */
	const ScheduleT* fBodySchedule; /**< body force schedule */
	dArrayT fBody; /**< body force vector */
	
	/** shape functions */
	MeshFreeNodalShapeFunctionT* fNodalShapes;
	
	/** underlying Element connectivities. Needed only for MLS stuff right now */
	ArrayT<const iArray2DT*> fElementConnectivities;

	/** equation numbers */
	RaggedArray2DT<int> fEqnos;

	/* traction data */
	dArray2DT fTractionVectors;
	iArrayT fTractionBoundaryCondition;

};

} /* namespace Tahoe */

#endif /* _SCNIMF_T_H_ */


