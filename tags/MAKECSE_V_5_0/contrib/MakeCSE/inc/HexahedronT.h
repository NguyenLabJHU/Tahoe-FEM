/*
 * File: HexahedronT.h
 */

/*
 * created      : PAK (10/22/97)
 * last modified: PAK (05/06/99)
 */

#ifndef _HEXADHEDRON_T_H_
#define _HEXADHEDRON_T_H_

/* base class */
#include "GeometryBaseT.h"

class HexahedronT: public GeometryBaseT 
{
  public:

	/* constructor */
	HexahedronT(int numnodes);

 	/* compute local shape functions and derivatives */
 	virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
 		dArrayT& weights);

   	/* set the values of the nodal extrapolation matrix */
 	virtual void SetExtrapolation(dMatrixT& extrap);

 	/* return the local node numbers for each facet of the element 
 	 * numbered to produce at outward normal in the order: vertex
 	 * nodes, mid-edge nodes, mid-face nodes */
	virtual void NodesOnFacet(int facet, iArrayT& facetnodes) const;
	virtual void NumNodesOnFacets(iArrayT& num_nodes) const;

	/* returns the nodes on each facet needed to determine neighbors
	 * across facets */
	virtual void NeighborNodeMap(iArray2DT& facetnodes) const;
	
	/* return geometry and number of nodes on each facet */
	virtual void FacetGeometry(ArrayT<GeometryCode>& facet_geom, 
		iArrayT& facet_nodes) const;
};

#endif /* _HEXADHEDRON_T_H_ */
