/* $Id: iGridManager2DT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/09/1997)                                          */
/* iNodeT grid                                                            */

#include "iGridManager2DT.h"
#include "iArrayT.h"

/* constructor */
iGridManager2DT::iGridManager2DT(int nx, int ny, const dArray2DT& coords,
	const iArrayT* nodes_used):
	GridManager2DT<iNodeT>(nx, ny, coords, nodes_used),
	fCoords(coords),
	fNodesUsed(nodes_used)
{

}	

/* neighbors - returns neighbors coords(n) (SELF not included) */
void iGridManager2DT::Neighbors(int n, double tol, AutoArrayT<int>& neighbors)
{
	/* initialize */
	neighbors.Allocate(0);
	
	/* fetch prospective neighbors */
	double* target = fCoords(n);
	const AutoArrayT<iNodeT>& hits =  HitsInRegion(target, tol);

	/* search through list */
	double tolsqr = tol*tol;
	int   thistag = n;
	for (int i = 0; i < hits.Length(); i++)
		if (hits[i].Tag() != thistag)
		{
			double* coords = hits[i].Coords();
			
			double dx = target[0] - coords[0];
			double dy = target[1] - coords[1];
			double dsqr = dx*dx + dy*dy;
			
			/* add to neighbor list */
			if (dsqr <= tolsqr) neighbors.Append(hits[i].Tag());
		}
}

/* reconfigure grid with stored coordinate data */
void iGridManager2DT::Reset(void)
{
	/* inherited */
	GridManager2DT<iNodeT>::Reset(fCoords, fNodesUsed);

	/* re-insert coordinate data */
	if (fNodesUsed == NULL)
	{	
		iNodeT temp;
		for (int i = 0; i < fCoords.MajorDim(); i++)
		{
			temp.Set(fCoords(i), i);
			Add(temp);
		}
	}
	else
	{
		iNodeT temp;
		int* dex = fNodesUsed->Pointer();
		for (int i = 0; i < fNodesUsed->Length(); i++)
		{
			int node = *dex++;
			temp.Set(fCoords(node), node);
			Add(temp);
		}
	}	
}
