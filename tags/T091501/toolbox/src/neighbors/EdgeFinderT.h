/* $Id: EdgeFinderT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (02/14/1998)                                          */
/* Class to determine element neighbors based on the connectivies.        */

#ifndef _EDGE_FINDER_T_H_
#define _EDGE_FINDER_T_H_

/* direct members */
#include "iArray2DT.h"
#include "iArrayT.h"
#include "RaggedArray2DT.h"

class EdgeFinderT
{
public:

	/* constructor */
	EdgeFinderT(const iArray2DT& connects, const iArray2DT& nodefacetmap);

	/* clear (and free) all data */
	void Clear(void);

	/* element edge data */
	const iArray2DT& Neighbors(void);	
	
	// Other additions?
	// (1) return list border nodes
	// (2) return specs {elem, facet} for border facets
	// (3) return specs {elem, facet} for internal facets

private:

	/* get dimensions from the connectivity set */
	void SetDimensions(void);

	/* set elements(node) data */
	void SetInverseConnects(void);

	/* find facet of elem_j that matches facet i of elem_i */
	int FindMatchingFacet(int facet, const int* elem_i,
		const int* elem_j) const;
		
protected:

	/* connectivities */
	const iArray2DT& fConnects;
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

#endif /* _EDGE_FINDER_T_H_ */
