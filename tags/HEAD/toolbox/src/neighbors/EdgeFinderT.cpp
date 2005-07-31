/* $Id: EdgeFinderT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (02/14/1998)                                          */
/* Class to determine element neighbors based on the connectivies.        */
/* The neighboring element numbers (taken from position in the list       */
/* of connectivities) has the same dimension as the number of element     */
/* nodes. The order of neighbors corresponds to the order of edges        */
/* as they appear in the connectivities. If no neighbor is found for      */
/* any edge, the neighbor number is -1.                                   */

#include "EdgeFinderT.h"

#include "AutoArrayT.h"

/* indices for data flags */
int kNumFlags   = 3;
int kDims       = 0;
int kNeighbor   = 1;
int kInvConnect = 2;

/* constructor */
EdgeFinderT::EdgeFinderT(const iArray2DT& connects,
	const iArray2DT& nodefacetmap):
	fCurrent(kNumFlags),
	fConnects(connects),
	fNumFacets(nodefacetmap.MajorDim()),
	fKeyNodes(nodefacetmap.Max() + 1), // assuming key nodes appear sequentially as
	                                   // the first entries in the connectivities
	fNodeFacetMap(nodefacetmap)
{
	/* check */
	if (fNumFacets < 3) throw eGeneralFail; // at least tri's?
	if (fKeyNodes < fNumFacets) throw eGeneralFail;

	/* initialize */
	Clear();
}

/* clear (and free) all data */
void EdgeFinderT::Clear(void)
{
	/* flag all not current */
	fCurrent = 0;

	/* set node number range */
	fMinNum   =-1;
	fMaxNum   =-1;
	fNumNodes = 0;
	
	/* release memory */
	fNeighbors.Free();
	fInvConnects.Free();
}

/* element edge data */
const iArray2DT& EdgeFinderT::Neighbors(void)
{
	if (!fCurrent[kNeighbor])
	{
		fCurrent[kNeighbor] = 1;

		/* get connectivity set dimensions */
		SetDimensions();

		/* set elements(node) */
		SetInverseConnects();
		
		/* allocate and initialize neighbor data */
		int nel = fConnects.MajorDim();
		fNeighbors.Allocate(nel, fNumFacets);
		fNeighbors = -1;

		/* work space */
		iArrayT hit_count(nel);
		hit_count = 0;
		
		/* set neighbors */
		int nfn = fNodeFacetMap.MinorDim();
		AutoArrayT<int> hit_elems;
		for (int i = 0; i < nel; i++)
		{
			int* neigh_i = fNeighbors(i);
			int* elem_i  = fConnects(i);
			for (int j = 0; j < fNumFacets; j++)
			{
				/* neighbor not set */
				if (*neigh_i == -1)
				{
					/* tag all neighboring elements */
					int* nodemap = fNodeFacetMap(j);
					for (int k = 0; k < nfn; k++)
					{
						/* location of elems(node) data */
						int row = elem_i[nodemap[k]] - fMinNum;
						
						/* the data */
						int* elems = fInvConnects(row);
						int    dim = fInvConnects.MinorDim(row);
						
						/* tag */
						for (int l = 0; l < dim; l++)
						{
							hit_elems.AppendUnique(*elems);
							hit_count[*elems++]++;
						}
					}
					
					/* find facet neighbor */
					int num_hit = hit_elems.Length();
					int neighbor = -1;
					for (int l = 0; l < num_hit && neighbor < 0; l++)
					{
						int elem_l = hit_elems[l];
						if (elem_l != i && hit_count[elem_l] == nfn)
							neighbor = elem_l;
					}
					
					/* found neighbor */
					if (neighbor != -1)
					{
						/* determine neighbor's facet */
						int facet_l = FindMatchingFacet(j, elem_i, fConnects(neighbor));
				
					#if __option(extended_errorcheck)
						/* neighbor's neighbor should be unset */
						if (fNeighbors(neighbor, facet_l) != -1)
						{
							cout << "\n EdgeFinderT::Neighbors: neighbor's neighbor already set" << endl;
							throw eGeneralFail;
						}
					#endif
						
						/* cross link */
						*neigh_i = neighbor;
						fNeighbors(neighbor, facet_l) = i;
					}
				}
				
				/* reset hit counts */
				int*   hits = hit_elems.Pointer();
				int num_hit = hit_elems.Length();
				for (int k = 0; k < num_hit; k++)
					hit_count[*hits++] = 0;
				
				/* initialize hit list */	
				hit_elems.Allocate(0);
		
				/* next facet */
				neigh_i++;
			}
		}
	}

	return fNeighbors;
}

/**********************************************************************
* Private
**********************************************************************/

/* get dimensions from the connectivity set */
void EdgeFinderT::SetDimensions(void)
{
	if (!fCurrent[kDims])
	{
		/* mark as current */
		fCurrent[kDims] = 1;
	
		/* check */
		if (fKeyNodes > fConnects.MinorDim()) throw eOutOfRange;
	
		/* set node number range */
		fMinNum   = fConnects.Min();
		fMaxNum   = fConnects.Max();
		fNumNodes = fMaxNum - fMinNum + 1;
	}
}

/* set elements(node) data */
void EdgeFinderT::SetInverseConnects(void)
{
	if (!fCurrent[kInvConnect])
	{
		/* mark as current */
		fCurrent[kInvConnect] = 1;
		
		/* workspace */
		AutoFill2DT<int> invconnects(fNumNodes, 25, fNumFacets);

		/* generate map */
		int  nen = fConnects.MinorDim();
		int* pel = fConnects(0);
		for (int i = 0; i < fConnects.MajorDim(); i++)
		{
			int* pel_i = pel;
			for (int j = 0; j < fKeyNodes; j++)
			{
				int row = (*pel_i++) - fMinNum;
				invconnects.AppendUnique(row, i);
			}

			pel += nen;
		}		

		/* copy data */
		fInvConnects.Copy(invconnects);
	}
}

/* find facet of elem_j that matches facet of elem_i */
int EdgeFinderT::FindMatchingFacet(int facet_i, const int* elem_i,
	const int* elem_j) const
{
	int nfn = fNodeFacetMap.MinorDim();
	int facet_j = -1;
	for (int j = 0; j < fNumFacets && facet_j < 0; j++)
	{
		int* nodemap_i = fNodeFacetMap(facet_i) + (nfn - 1);
		int* nodemap_j = fNodeFacetMap(j);

		/* find starting point */
		int node_i = elem_i[*nodemap_i--];
		int ji = -1;
		for (int k = 0; k < nfn && ji < 0; k++)
			if (elem_j[*nodemap_j++] == node_i) ji = k;
	
		/* found shared node */
		if (ji > -1)
		{
			int match = 1;
			
			/* check match of remaining facet nodes */
			for (int k = 1; k < nfn && match == 1; k++)
			{
				/* wrap list */
				if (++ji == nfn) nodemap_j -= nfn;
			
				/* compare */
				if (elem_j[*nodemap_j++] != elem_i[*nodemap_i--]) match = 0;
			}
			
			if (match == 1) facet_j = j;
		}
	}

#if __option(extended_errorcheck)
	/* facet not found */
	if (facet_j == -1)
	{
		cout << "\n EdgeFinderT::FindMatchingFacet: failed" << endl;
		throw eGeneralFail;
	}
#endif

	return facet_j;	
}
