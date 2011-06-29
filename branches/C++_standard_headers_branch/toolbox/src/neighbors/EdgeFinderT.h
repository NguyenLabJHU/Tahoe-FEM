/* $Id: EdgeFinderT.h,v 1.6 2003-11-21 22:41:59 paklein Exp $ */
/* created: paklein (02/14/1998) */
#ifndef _EDGE_FINDER_T_H_
#define _EDGE_FINDER_T_H_

/* direct members */
#include "iArray2DT.h"
#include "iArrayT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

/** class to determine element neighbors based on the connectivies.
 * The neighboring element numbers (taken from position in the list
 * of connectivities) has the same dimension as the number of element
 * nodes. The order of neighbors corresponds to the order of edge
 * as they appear in the connectivities. If no neighbor is found for
 * any edge, the neighbor number is -1. 
 * \note Currently, all element blocks passed in must have the same
 *       topology, i.e., number of element nodes. */
class EdgeFinderT
{
public:

	/** constructor
	 * \param connects list of connectivity blocks to consider as a single body
	 * \param nodefacetmap list of the local node numbers on each face of the
	 *        element*/
	EdgeFinderT(const ArrayT<const iArray2DT*>& connects, const iArray2DT& nodefacetmap);

	/** clear (and free) all data */
	void Clear(void);

	/** determine element edge data. Elements numbered sequentially through all blocks. */
	const iArray2DT& Neighbors(void);	
	
	/** return the "bounding" elements. For the body comprised of the given
	 * list of element blocks, determine the elements with faces on the body
	 * boundary and the corresponding neighbors. All returned arrays are
	 * dimensioned during the call.
	 * \param elements elements within each block that have faces on the body 
	 *        boundary. Element are numbered sequentially through all blocks
	 * \param neighbors neighoring elements for bounding element 
	 */
	void BoundingElements(iArrayT& elements, iArray2DT& neighbors);

	/** element faces on the group "surface".
	 * \param surface_facets element faces defining the body boundary: [nf] x [nfn]
	 * \param surface_nodes nodes comprising the surface facets
	 */
	void SurfaceFacets(iArray2DT& surface_facets, iArrayT& surface_nodes);
	void SurfaceFacets(iArray2DT& surface_facets, iArrayT& surface_nodes, 
					   iArrayT& facet_numbers, iArrayT& elem_numbers);

	/** element faces on the group "surface" grouped into contiguous patches.
	 * Determines the surfaces faces and groups them into contiguous patches based 
	 * on the connectivities. Also, returns the nodes comprising the facets.
	 */
	void SurfaceFacets(ArrayT<iArray2DT>& surface_facet_sets, 
		iArrayT& surface_nodes);
		
	/** return the nodes lying on the surface of the body */
	void SurfaceNodes(iArrayT& surface_nodes);
	
	/** map element index to pointer to element nodes across all connectivities
	 * being operated on */
	const int* ElementNodes (int index) const;
	
	// Other additions?
	// (1) return list border nodes
	// (2) return specs {elem, facet} for border facets
	// (3) return specs {elem, facet} for internal facets

private:

	/** get dimensions from the connectivity set */
	void SetDimensions(void);

	/** set elements(node) data */
	void SetInverseConnects(void);

	/** find facet of elem_j that matches facet i of elem_i */
	int FindMatchingFacet(int facet, const int* elem_i,
		const int* elem_j) const;
		
protected:

	/* connectivities */
	ArrayT<const iArray2DT*> fConnects;
	iArrayT fStartNumber;
	int fNumElements;
	int fNumFacets;
	int fKeyNodes;

	/* nodes(facet) map */
	const iArray2DT fNodeFacetMap;

	/* range of node numbers */
	int fMinNum;
	int fMaxNum;
	int fNumNodes;

	/* data flags */
	iArrayT fCurrent;
	
	/* data */
	iArray2DT fNeighbors;             // 0: element neighbor lists
	RaggedArray2DT<int> fInvConnects; // 1: elements(node)
};

} // namespace Tahoe 
#endif /* _EDGE_FINDER_T_H_ */
