#ifndef _VORONOI_DIAGRAM_T_
#define _VORONOI_DIAGRAM_T_

/* base class */
#include "CellGeometryT.h"

#include "LinkedListT.h"

#ifdef __QHULL__
#include "CompGeomT.h"
#endif

namespace Tahoe {

/** base class for particle types */
class VoronoiDiagramT: public CellGeometryT
{
public:

	/** constructor */
	VoronoiDiagramT(const ElementSupportT& support, bool isAxisymmetric);
	VoronoiDiagramT(void);

	/** destructor */
	~VoronoiDiagramT(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected: /* for derived classes only */

	/** \name types needed for the Voronoi diagram calculation */
	/*@{*/
#ifndef __QHULL__
	/** Basic structure -- hullMap[i][j] gives index of jth vertex of structure i */
	typedef ArrayT<iArrayT> ConvexHullMap;

	/** Voronoi diagram facet structure is an array of convex hulls  */
	typedef ArrayT< ConvexHullMap > VoronoiDiagramMap;
#endif
	/*@}*/

	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);
	
	/** transfers data from QHULL and computes new data structures. This function 
		initializes fNonDeloneEdges, fNonDeloneNormals, fSelfDualFacets, fBoundaryIntegrationWeights,
		and fVoronoiCellCentroids. If qhull is not used, these data structures are
		read in by VoronoiDiagramFromFile */
	void InitializeVoronoiData(void);
	
	/** generate data structures for integration over the body boundary */
	virtual void BoundaryShapeFunctions(RaggedArray2DT<double>& phis, RaggedArray2DT<int>& supports, dArray2DT& normals);
	
	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(RaggedArray2DT<int>& nodalCellSupports, RaggedArray2DT<dArrayT>& bVectorArray,
								  dArrayT& cellVolumes, dArray2DT& cellCentroids, RaggedArray2DT<double>& circumferential_B);
	
	/** write out Voronoi diagram data. This function is only called when qhull is used
		to compute the clipped Voronoi diagram. */
	void VoronoiDiagramToFile(ofstreamT& vout);
	
	/** read in Voronoi diagram data. This function reads in data structures when qhull
		is not used to create them. It initializes fVoronoiVertices, fCellVolumes,
		fDeloneEdges, fDualFacets, fNonDeloneEdges, fNonDeloneNormals, fSelfDualFacets, 
		fBoundaryIntegrationWeights, and fVoronoiCellCentroids. If qhull is used, these
		data structures are either gotten directly from qhull or created in 
		InitializeVoronoiData. */
	void VoronoiDiagramFromFile(ifstreamT& vin);
	
	/** For 3D only, identifies nodes at corners of the body and the edges
            of the body. Currently assumes that elements meet the body boundary
	    only by faces--i.e. there are no degenerate intersections between
	    elements and the element block boundary. */
	void FindCornersAndEdges(iArrayT& boundaryFacets, iArrayT& boundaryElements);

protected:
	
	/** \name Data Structures for Voronoi Decomposition */
	/*@{*/
	
	/** these are dual to Voronoi facets. They have minor dimension of 2 . 
	     Difference in the two points is parallel to the normal vector of the 
	     Voronoi facet dual to the Delone edge. This data structure is created
	     by qhull or read in from a text file. */
	iArray2DT fDeloneEdges;

	/** Voronoi facets dual to the Delone Edges -- This data structure is currently
		a list of vertices of the facets. In 2D, all facets are simplicial and this
		is fine. In 3D, it will change. */
	RaggedArray2DT<int> fDualFacets; 
	
	/** boundary facets. See qhull/CompGeomT.h for a definition of the self-dual
		terminology -- this data structure is 
		created in InitializeVoronoiData based on qhull's clipped Voronoi diagram. It is then
		written to the geometry text file. This is essentially a flattened version of
		qhulls selfDual data structure */
	iArray2DT fSelfDualFacets; 

	/** connectivity of boundary nodes. Currently determined from an underlying 
	    element connectivity */
	iArray2DT fBoundaryConnectivity; 
	
	/** union of nodes in fBoundaryConnectivity */
	iArrayT fBoundaryNodes; 
	
	/** true if boundary connectivity is simplicial */
	bool fBoundaryIsTriangulated; 
	
	/** additional edges associated only with one node -- this data structure is 
		created in InitializeVoronoiData based on qhull's clipped Voronoi diagram. It is then
		written to the geometry text file. Since these facets are self dual, there is only
		one Delone vertex per facet, and this vertex is the boundary vertex. The list of
		boundary vertices for each self-dual facet is contained in this array of integers. */
	iArrayT fNonDeloneEdges; 
	
	/** normal vectors of the facets for those edges -- this data structure is 
		created in InitializeVoronoiData based on qhull's clipped Voronoi diagram. It is then
		written to the geometry text file.*/
	dArray2DT fNonDeloneNormals;

	/** Compute or read the Voronoi Diagram -- these data structures are created by qhull */	
#ifdef __QHULL__	
	CompGeomT* fVoronoi;
#else
	void* fVoronoi;
#endif
	bool qComputeVoronoiCell, qJustVoronoiDiagram;
	StringT vCellFile;
	
	/** The coordinates of the Voronoi diagram after it has been clipped by the body boundary */
	dArray2DT fVoronoiVertices; 
	
	/** The volumes of the Voronoi Cells */
	dArrayT fVoronoiCellVolumes;
	
	/** The centroids of the Voronoi Cells */
	dArray2DT fVoronoiCellCentroids;
	
	/** */
	dArrayT fBoundaryIntegrationWeights;
	
	/* These should be local, but right now they're only here b/c of the kludgy way of the boundary handling */
	/** temporary storage for integration over the body boundary */
	ArrayT< LinkedListT<double> > boundary_phi;
	ArrayT< LinkedListT<int> > boundary_supports;
	
};

} /* namespace Tahoe */


#endif /* _VORONOI_DIAGRAM_T_ */

