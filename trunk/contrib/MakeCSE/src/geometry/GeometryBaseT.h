/*
 * File: GeometryBaseT.h
 *
 * Class to initialize shape function arrays with
 * geometry specific values.
 *
 */

/*
 * created      : PAK (10/21/1997)
 * last modified: PAK (10/10/1999)
 */

#ifndef _GEOMETRY_BASE_T_H_
#define _GEOMETRY_BASE_T_H_

/* base class */
#include "GeometryT.h"

/* forward declarations */
class dArray2DT;
template <class TYPE> class ArrayT;
class iArrayT;
class dArrayT;
class dMatrixT;
class iArray2DT;
template <class TYPE> class ArrayT;

class GeometryBaseT: public GeometryT
{
  public:

	/* constructor */
	GeometryBaseT(int numnodes, int numfacets);

	/* destructor */
	virtual ~GeometryBaseT(void);

	/* returns the number of element facets */
	int NumFacets(void) const;

 	/* compute local shape functions and derivatives */
 	virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
 		dArrayT& weights) = 0;

   	/* set the values of the nodal extrapolation matrix */
 	virtual void SetExtrapolation(dMatrixT& extrap) = 0;
 	
 	/* return the local node numbers for each facet of the element 
 	 * numbered to produce at outward normal in the order: vertex
 	 * nodes, mid-edge nodes, mid-face nodes */
	virtual void NodesOnFacet(int facet, iArrayT& facets) const = 0;
	virtual void NumNodesOnFacets(iArrayT& num_nodes) const = 0;
//	virtual void NodesOnFacet(RaggedArrayT<int>& facets) const = 0;

	/* returns the nodes on each facet needed to determine neighbors
	 * across facets */
	virtual void NeighborNodeMap(iArray2DT& facetnodes) const = 0;

	/* return geometry and number of nodes on each facet */
	virtual void FacetGeometry(ArrayT<GeometryCode>& facet_geom, 
		iArrayT& facet_nodes) const = 0;
		
  protected:
  
  	/* number of domain nodes */
  	int fNumNodes;
  	int fNumFacets;
};

/* returns the number of element facets */
inline int GeometryBaseT::NumFacets(void) const { return fNumFacets; }

#endif /* _GEOMETRY_BASE_T_H_ */
