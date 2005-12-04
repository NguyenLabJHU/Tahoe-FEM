/* $Id: PentahedronT.h,v 1.4 2005-12-04 16:56:29 paklein Exp $ */
/* created: sawimme (10/22/1999) */
#ifndef _PENTAHEDRON_T_H_
#define _PENTAHEDRON_T_H_

/* base class */
#include "GeometryBaseT.h"

namespace Tahoe {

/** 3D pentadedral parent domain. Normals on triangle facets point outward, 
 * others point inward */
class PentahedronT: public GeometryBaseT
{
public:

	/** constructor */
	PentahedronT(int numnodes);

	/** return the geometry code */
	virtual GeometryT::CodeT Geometry(void) const { return kPentahedron; };

	/** evaluate the shape functions. See 
	 * GeometryBaseT::EvaluateShapeFunctions for documentation */
	virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const;

	/** evaluate the shape functions and gradients. See 
	 * GeometryBaseT::EvaluateShapeFunctions for documentation */
	virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, 
		dArray2DT& DNa) const;

	/** evaluate the shape functions and gradients. See 
	 * GeometryBaseT::SetLocalShape for documentation */
	virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
		dArrayT& weights) const;

	/* set the values of the nodal extrapolation matrix */
	virtual void SetExtrapolation(dMatrixT& extrap) const;

	/* return the local node numbers for each facet of the element
	 * numbered to produce at outward normal in the order: vertex
	 * nodes, mid-edge nodes, mid-face nodes */
	virtual void NodesOnFacet(int facet, iArrayT& facetnodes) const;
	virtual void NumNodesOnFacets(iArrayT& num_nodes) const;

	/** return the local node numbers for each edge of element where
	 * the row number corresponds with the canonical numbering of the
	 * edges in the element. Edges are only defined for 3D domain geometries. */
	virtual void NodesOnEdges(iArray2DT& nodes_on_edges) const;

	/* returns the nodes on each facet needed to determine neighbors
	 * across facets */
	virtual void NeighborNodeMap(iArray2DT& facetnodes) const;
	
	/* return geometry and number of nodes on each facet */
	virtual void FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const;
};

} // namespace Tahoe 
#endif /* _PENTAHEDRON_T_H_ */
