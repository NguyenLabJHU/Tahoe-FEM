/* $Id: EdgeFinderT.cpp,v 1.9 2003-11-21 22:41:59 paklein Exp $ */
/* created: paklein (02/14/1998) */
#include "EdgeFinderT.h"
#include "AutoArrayT.h"
#include "GraphT.h"

using namespace Tahoe;

/* indices for data flags */
int kNumFlags   = 3;
int kDims       = 0;
int kNeighbor   = 1;
int kInvConnect = 2;

/* constructor */
EdgeFinderT::EdgeFinderT(const ArrayT<const iArray2DT*>& connects,
	const iArray2DT& nodefacetmap):
	fConnects(connects.Length()),
	fStartNumber (connects.Length()),
	fNumElements(0),
	fNumFacets(nodefacetmap.MajorDim()),
	fKeyNodes(nodefacetmap.Max() + 1), // assuming key nodes appear sequentially as
	                                   // the first entries in the connectivities
	fNodeFacetMap(nodefacetmap),
	fCurrent(kNumFlags)
{
	/* check */
	if (fNumFacets < 3) throw ExceptionT::kGeneralFail; // at least tri's?
	if (fKeyNodes < fNumFacets) throw ExceptionT::kGeneralFail;

	for (int i=0; i < connects.Length(); i++)
	  {
	    fStartNumber [i] = fNumElements;
	    fConnects[i] = connects[i];
	    fNumElements += connects[i]->MajorDim();
	  }

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
		fNeighbors.Dimension(fNumElements, fNumFacets);
		fNeighbors = -1;

		/* work space */
		iArrayT hit_count(fNumElements);
		hit_count = 0;
		
		/* set neighbors */
		int nfn = fNodeFacetMap.MinorDim();
		AutoArrayT<int> hit_elems;
		for (int i = 0; i < fNumElements; i++)
		{
			int* neigh_i = fNeighbors(i);
			const int* elem_i  = ElementNodes(i);
			for (int j = 0; j < fNumFacets; j++)
			{
				/* neighbor not set */
				if (*neigh_i == -1)
				{
					/* tag all neighboring elements */
					const int* nodemap = fNodeFacetMap(j);
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
						int facet_l = FindMatchingFacet(j, elem_i, ElementNodes(neighbor));
				
					#if __option(extended_errorcheck)
						/* neighbor's neighbor should be unset */
						if (fNeighbors(neighbor, facet_l) != -1)
						{
							cout << "\n EdgeFinderT::Neighbors: neighbor's neighbor already set" << endl;
							throw ExceptionT::kGeneralFail;
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
				hit_elems.Dimension(0);
		
				/* next facet */
				neigh_i++;
			}
		}
	}

	return fNeighbors;
}

/* return the "bounding" elements */
void EdgeFinderT::BoundingElements(iArrayT& elements, iArray2DT& neighbors)
{
	/* element neighbor list */
	const iArray2DT& all_neighbors = Neighbors();

	/* determine total number of elements */
	int nel = 0;
	for (int i = 0; i < fConnects.Length(); i++)
		nel += fConnects[i]->MajorDim();

	/* collect list of bounding elements */
	AutoArrayT<int> borders;
	iArrayT element;
	for (int i = 0; i < nel; i++)
	{
		all_neighbors.RowAlias(i, element);
	
		/* has "free" edge */
		if (element.HasValue(-1)) borders.Append(i);
	}
	elements.Dimension(borders.Length());
	borders.CopyInto(elements);
	
	/* copy bounding element neighbor lists */
	neighbors.Dimension(elements.Length(), all_neighbors.MinorDim());
	neighbors.RowCollect(elements, all_neighbors);
}

/* surface facets */
void EdgeFinderT::SurfaceFacets(iArray2DT& surface_facets, iArrayT& surface_nodes)
{
	/* find bounding elements */
	iArrayT   border_elems;
	iArray2DT border_neighs;
	BoundingElements(border_elems, border_neighs);
		
	/* collect nodes on facets info */
//TEMP
#if 0
	ArrayT<iArrayT> facetnodemap(geometry->NumFacets());
	for (int i2 = 0; i2 < facetnodemap.Length(); i2++)
		geometry->NodesOnFacet(i2, facetnodemap[i2]);	
#endif

	/* collect surface facets (with "outward" normal ordering) */
	AutoArrayT<int> border_nodes;
	int surf_count = 0;
	int num_facets = fNodeFacetMap.MajorDim();
	int num_facet_nodes = fNodeFacetMap.MinorDim();
	border_nodes.Dimension(0);
	surface_facets.Dimension(border_neighs.Count(-1), num_facet_nodes);
	for (int i = 0; i < border_elems.Length(); i++)
	{
		/* element connectivity */
		const int* elem = ElementNodes(border_elems[i]);

		/* find open sides */
		int found_open = 0;
		int* pneigh = border_neighs(i);
		for (int j = 0; j < num_facets; j++)
		{
			/* open face */
			if (*pneigh == -1)
			{
				/* set flag */
				found_open = 1;
				
				/* collect facet nodes */
				int* pfacet = surface_facets(surf_count++);
				const int* facet_nodes = fNodeFacetMap(j);
				for (int k = 0; k < num_facet_nodes; k++)
				{
					int node = elem[*facet_nodes++];
					*pfacet++ = node;
					border_nodes.AppendUnique(node);
					// better just to keep a "nodes used" map?
				}
			}	
			pneigh++;
		}
	
		/* no open facet */	
		if (!found_open)
		{
			cout << "\n EdgeFinderT::SurfaceFacets: error building surface facet list" << endl;
			throw ExceptionT::kGeneralFail;
		}	
	}

	/* return value */
	surface_nodes.Dimension(border_nodes.Length());
	border_nodes.CopyInto(surface_nodes);
}
void EdgeFinderT::SurfaceFacets(iArray2DT& surface_facets, 
	iArrayT& surface_nodes, iArrayT& facet_numbers, iArrayT& elem_numbers)
{
	/* find bounding elements */
	iArrayT   border_elems;
	iArray2DT border_neighs;
	BoundingElements(border_elems, border_neighs);
		
	/* collect nodes on facets info */
//TEMP
#if 0
	ArrayT<iArrayT> facetnodemap(geometry->NumFacets());
	for (int i2 = 0; i2 < facetnodemap.Length(); i2++)
		geometry->NodesOnFacet(i2, facetnodemap[i2]);	
#endif

	/* collect surface facets (with "outward" normal ordering) */
	AutoArrayT<int> border_nodes;
	int surf_count = 0;
	int num_facets = fNodeFacetMap.MajorDim();
	int num_facet_nodes = fNodeFacetMap.MinorDim();
	border_nodes.Dimension(0);
	surface_facets.Dimension(border_neighs.Count(-1), num_facet_nodes);
	facet_numbers.Dimension(surface_facets.MajorDim());
	elem_numbers.Dimension(surface_facets.MajorDim());
	for (int i = 0; i < border_elems.Length(); i++)
	{
		/* element connectivity */
		const int* elem = ElementNodes(border_elems[i]);

		/* find open sides */
		int found_open = 0;
		int* pneigh = border_neighs(i);
		for (int j = 0; j < num_facets; j++)
		{
			/* open face */
			if (*pneigh == -1)
			{
				/* set flag */
				found_open = 1;
				
				/* collect facet nodes */
				const int* facet_nodes = fNodeFacetMap(j);
				facet_numbers[surf_count] = j; 			
				elem_numbers[surf_count] = border_elems[i];	
				int* pfacet = surface_facets(surf_count++);
				for (int k = 0; k < num_facet_nodes; k++)
				{
					int node = elem[*facet_nodes++];
					*pfacet++ = node;
					border_nodes.AppendUnique(node);
					// better just to keep a "nodes used" map?
				}
			}	
			pneigh++;
		}
	
		/* no open facet */	
		if (!found_open)
		{
			cout << "\n EdgeFinderT::SurfaceFacets: error building surface facet list" << endl;
			throw ExceptionT::kGeneralFail;
		}	
	}

	/* return value */
	surface_nodes.Dimension(border_nodes.Length());
	border_nodes.CopyInto(surface_nodes);
}
/* with surface facets sorted into connected sets */
void EdgeFinderT::SurfaceFacets(ArrayT<iArray2DT>& surface_facet_sets,
	iArrayT& surface_nodes)
{
	/* collect all surface facets */
	iArray2DT surface_facets;
	SurfaceFacets(surface_facets, surface_nodes);

	/* graph object */
	GraphT graph;
	graph.AddGroup(surface_facets);
	graph.MakeGraph();

	iArrayT branch_map;
	graph.LabelBranches(surface_nodes, branch_map);
	
	/* sort surfaces */
	int num_branches = branch_map.Max() + 1;
	surface_facet_sets.Dimension(num_branches);
	if (num_branches == 1)
		surface_facet_sets[0] = surface_facets;
	else
	{
		/* surfaces in each set */
		iArrayT count(num_branches);
		int size = surface_facets.MinorDim();
		count = 0;
		int* psurf = surface_facets(0);
		for (int i = 0; i < surface_facets.MajorDim(); i++)
		{
			count[branch_map[*psurf]]++;
			psurf += size;
		}

		/* set to surfaces map */
		RaggedArray2DT<int> set_data;
		set_data.Configure(count);

		count = 0;
		psurf = surface_facets(0);
		for (int j = 0; j < surface_facets.MajorDim(); j++)
		{
			int branch = branch_map[*psurf];
			*(set_data(branch) + count[branch]) = j;

			count[branch]++;
			psurf += size;
		}
			
		/* copy in */
		for (int k = 0; k < num_branches; k++)
		{
			surface_facet_sets[k].Dimension(set_data.MinorDim(k), size);
			surface_facet_sets[k].RowCollect(set_data(k), surface_facets);
		}
	}
}		

/* with surface facets sorted into connected sets */
void EdgeFinderT::SurfaceNodes(iArrayT& surface_nodes)
{
	/* collect all surface facets */
	iArray2DT surface_facets;
	SurfaceFacets(surface_facets, surface_nodes);
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
		int nen = fConnects[0]->MinorDim();
		for (int i=0; i < fConnects.Length(); i++)
		  {
		    if (fKeyNodes > fConnects[i]->MinorDim()) throw ExceptionT::kOutOfRange;
		    if (nen != fConnects[i]->MinorDim()) throw ExceptionT::kSizeMismatch;
		  }

		/* set node number range */
		iArrayT mins (fConnects.Length());
		iArrayT maxes (fConnects.Length());
		for (int i=0; i < fConnects.Length(); i++)
		  {
		    mins[i] = fConnects[i]->Min();
		    maxes[i] = fConnects[i]->Max();
		  }

		fMinNum   = mins.Min();
		fMaxNum   = maxes.Max();
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
		AutoFill2DT<int> invconnects(fNumNodes, 1, 25, fNumFacets);

		/* generate map */
//		int  nen = fConnects[0]->MinorDim();
		for (int i = 0; i < fNumElements; i++)
		{
			const int* pel = ElementNodes(i);
			for (int j = 0; j < fKeyNodes; j++)
			{
				int row = (*pel++) - fMinNum;
				invconnects.AppendUnique(row, i);
			}
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
		const int* nodemap_i = fNodeFacetMap(facet_i) + (nfn - 1);
		const int* nodemap_j = fNodeFacetMap(j);

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
		throw ExceptionT::kGeneralFail;
	}
#endif

	return facet_j;	
}

const int* EdgeFinderT::ElementNodes (int index) const
{
  if (index < 0 || index >= fNumElements) throw ExceptionT::kOutOfRange;

  /* find the block */
  int block = 0;
//  int offset = 0;
  while (block+1 < fStartNumber.Length() && index >= fStartNumber[block+1])
    {
      block++;
      if (block > fConnects.Length()) throw ExceptionT::kOutOfRange;
    }

  int localindex = index - fStartNumber[block];
  const iArray2DT* conn = fConnects[block];
  return (*conn)(localindex);
}
