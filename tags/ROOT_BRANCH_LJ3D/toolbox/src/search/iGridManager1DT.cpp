/* $Id: iGridManager1DT.cpp,v 1.5 2003-11-21 22:42:07 paklein Exp $ */
#include "iGridManager1DT.h"
#include "iArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
iGridManager1DT::iGridManager1DT(int nx, const dArray2DT& coords,
	const ArrayT<int>* nodes_used):
	GridManager1DT<iNodeT>(nx, coords, nodes_used),
	fCoords(coords),
	fNodesUsed(nodes_used)
{

}	

/* neighbors - returns neighbors coords(n) (SELF not included) */
void iGridManager1DT::Neighbors(int n, double tol, AutoArrayT<int>& neighbors)
{
	/* initialize */
	neighbors.Dimension(0);
	
	/* fetch prospective neighbors */
	const double* target = fCoords(n);
	const AutoArrayT<iNodeT>& hits =  HitsInRegion(target, tol);

	/* search through list */
	double tolsqr = tol*tol;
	int   thistag = n;
	for (int i = 0; i < hits.Length(); i++)
		if (hits[i].Tag() != thistag)
		{
			const double* coords = hits[i].Coords();
			
			double dx = target[0] - coords[0];
			double dsqr = dx*dx;
			
			/* add to neighbor list */
			if (dsqr <= tolsqr) neighbors.Append(hits[i].Tag());
		}
}

void iGridManager1DT::Neighbors(int n, const ArrayT<double>& tol_xy, 
	AutoArrayT<int>& neighbors)
{
	/* check */
	if (tol_xy.Length() != 1)
	{
		cout << "\n iGridManager1DT::Neighbors: expecting tolerance list length 1: " 
		     << tol_xy.Length() << endl;
		throw ExceptionT::kSizeMismatch;
	}

	/* initialize */
	neighbors.Dimension(0);
	
	/* fetch prospective neighbors */
	const double* target = fCoords(n);
	const AutoArrayT<iNodeT>& hits =  HitsInRegion(target, tol_xy);

	/* search through list */
	double tol_x = tol_xy[0];
	int   thistag = n;
	for (int i = 0; i < hits.Length(); i++)
		if (hits[i].Tag() != thistag)
		{
			const double* coords = hits[i].Coords();
			
			double dx = fabs(target[0] - coords[0]);
			
			/* add to neighbor list */
			if (dx <= tol_x) neighbors.Append(hits[i].Tag());
		}
}

/* reconfigure grid with stored coordinate data */
void iGridManager1DT::Reset(void)
{
	/* inherited */
	GridManager1DT<iNodeT>::Reset(fCoords, fNodesUsed);

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
		const int* dex = fNodesUsed->Pointer();
		for (int i = 0; i < fNodesUsed->Length(); i++)
		{
			int node = *dex++;
			temp.Set(fCoords(node), node);
			Add(temp);
		}
	}	
}
