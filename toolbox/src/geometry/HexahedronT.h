/* $Id: HexahedronT.h,v 1.8 2005-12-04 16:56:29 paklein Exp $ */
/* created: paklein (10/22/1997) */
#ifndef _HEXADHEDRON_T_H_
#define _HEXADHEDRON_T_H_

/* base class */
#include "GeometryBaseT.h"

namespace Tahoe {

/** calculations over a hexahedral parent domain. */
class HexahedronT: public GeometryBaseT
{
public:

	/** constructor */
	HexahedronT(int numnodes);

	/** return the geometry code */
	virtual GeometryT::CodeT Geometry(void) const { return kHexahedron; };

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

	/* compute gradients of the "bubble" modes */
	virtual void BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const;

	/* set the values of the nodal extrapolation matrix */
	virtual void SetExtrapolation(dMatrixT& extrap) const;

	/** integration point gradient matrix */
	virtual void IPGradientTransform(int ip, dMatrixT& transform) const;

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
	virtual void FacetGeometry(ArrayT<CodeT>& facet_geom,
		iArrayT& facet_nodes) const;

	/** return true if the given point is within the domain defined by
	 * the list of coordinates.
	 * \param coords list of coordinates defining the domain
	 * \param point test point coordinates */
	virtual bool PointInDomain(const LocalArrayT& coords, const dArrayT& point) const;

	/** \name nodal subdomains, see GeometryBaseT for more information */
	/*@{*/
	/** subdomain geometry */
	virtual GeometryT::CodeT NodalSubDomainGeometry(void) const;

	/** number of nodes defining the nodal subdomain */
	virtual int NodalSubDomainNumPoints(void) const;
	
	/** compute the coordinates of the points defining the nodal subdomain */
	virtual void NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
		LocalArrayT& subdomain_coords) const;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _HEXADHEDRON_T_H_ */
