/* $Id: iGridManager3DT.cpp,v 1.6 2003-11-21 22:42:07 paklein Exp $ */
/* created: paklein (12/09/1997) */
#include "iGridManager3DT.h"
#include "iArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
iGridManager3DT::iGridManager3DT(int nx, int ny, int nz,
	const dArray2DT& coords, const ArrayT<int>* nodes_used):
	GridManager3DT<iNodeT>(nx, ny, nz, coords, nodes_used),
	fCoords(coords),
	fNodesUsed(nodes_used)
{

}	

/* neighbors - returns neighbors coords(n) (SELF not included) */
void iGridManager3DT::Neighbors(int n, double tol, AutoArrayT<int>& neighbors)
{
	/* initialize */
	neighbors.Dimension(0);
	
	/* fetch prospective neighbors */
	const double* target = fCoords(n);
	const AutoArrayT<iNodeT>& hits =  HitsInRegion(target,tol);

	/* search through list */
	double tolsqr = tol*tol;
	int   thistag = n;
	for (int i = 0; i < hits.Length(); i++)
		if (hits[i].Tag() != thistag)
		{
			const double* coords = hits[i].Coords();
			
			double dx = target[0] - coords[0];
			double dy = target[1] - coords[1];
			double dz = target[2] - coords[2];
			double dsqr = dx*dx + dy*dy + dz*dz;
			
			/* add to neighbor list */
			if (dsqr <= tolsqr) neighbors.Append(hits[i].Tag());
		}
}

void iGridManager3DT::Neighbors(int n, const ArrayT<double>& tol_xyz, AutoArrayT<int>& neighbors)
{
	/* initialize */
	neighbors.Dimension(0);
	
	/* fetch prospective neighbors */
	const double* target = fCoords(n);
	const AutoArrayT<iNodeT>& hits =  HitsInRegion(target, tol_xyz);

	/* search through list */
	double tol_x = tol_xyz[0];
	double tol_y = tol_xyz[1];
	double tol_z = tol_xyz[2];
	int   thistag = n;
	for (int i = 0; i < hits.Length(); i++)
		if (hits[i].Tag() != thistag)
		{
			const double* coords = hits[i].Coords();
			
			double dx = fabs(target[0] - coords[0]);
			double dy = fabs(target[1] - coords[1]);
			double dz = fabs(target[2] - coords[2]);
			
			/* add to neighbor list */
			if (dx <= tol_x && 
			    dy <= tol_y &&
			    dz <= tol_z) neighbors.Append(hits[i].Tag());
		}
}

/* reconfigure grid with stored coordinate data */
void iGridManager3DT::Reset(void)
{
	/* inherited */
	GridManager3DT<iNodeT>::Reset(fCoords, fNodesUsed);

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

/**********************************************************************
* Protected
**********************************************************************/

/* process hit list - returns all, and repeated, tags != skiptag */
void iGridManager3DT::ProcessHits(double* target, double tol, int skiptag,
	int& count, iArrayT& neighbors)
{
	const AutoArrayT<iNodeT>& hits = HitsInRegion(target, tol);

	/* search through list */
	double tolsqr = tol*tol;
	for (int i = 0; i < hits.Length(); i++)
		if (hits[i].Tag() != skiptag)
		{
			const double* coords = hits[i].Coords();
			
			double dx = target[0] - coords[0];
			double dy = target[1] - coords[1];
			double dz = target[2] - coords[2];
			double dsqr = dx*dx + dy*dy + dz*dz;
			
			/* add to neighbor list */
			if (dsqr <= tolsqr)
			{
				/* quit if list full */
				if (count == neighbors.Length()) throw ExceptionT::kGeneralFail;

				/* append to list */
				neighbors[count++] = hits[i].Tag();			
			}
		}
}	

/* process hit list - only returns unique tags > skiptag */
void iGridManager3DT::ProcessHitsSorted(double* target, double tol,
	int skiptag, AutoArrayT<int>& neighbors)
{
	const AutoArrayT<iNodeT>& hits = HitsInRegion(target, tol);

	/* search through list */
	double tolsqr = tol*tol;
	for (int i = 0; i < hits.Length(); i++)
	{
		/* current tag */
		int tag = hits[i].Tag();
	
		if (tag > skiptag)
		{
			const double* coords = hits[i].Coords();
			
			double dx = target[0] - coords[0];
			double dy = target[1] - coords[1];
			double dz = target[2] - coords[2];
			double dsqr = dx*dx + dy*dy + dz*dz;
			
			/* add to neighbor list */
			if (dsqr <= tolsqr) neighbors.AppendUnique(tag);
		}
	}
}	
