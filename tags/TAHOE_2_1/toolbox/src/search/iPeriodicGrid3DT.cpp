/* $Id: iPeriodicGrid3DT.cpp,v 1.6 2003-11-21 22:42:07 paklein Exp $ */
/* created: paklein (12/18/1997)                                          */

#include "iPeriodicGrid3DT.h"
#include "dArrayT.h"
#include "iArrayT.h"

/* constructor */

using namespace Tahoe;

iPeriodicGrid3DT::iPeriodicGrid3DT(int nx, int ny, int nz,
	const dArray2DT& coords, const ArrayT<int>* nodes_used,
	const dArrayT& periodicity):
	iGridManager3DT(nx, ny, nz, coords, nodes_used),
	fPeriodicity(periodicity)
{

}	

/* neighbors - returns the number of neighbors */
int iPeriodicGrid3DT::PeriodicNeighbors(int n, double tol, iArrayT& neighbors)
{
	/* initialize */
	fSortedHits.Dimension(0);
	int skiptag = n;

	double h_x = fPeriodicity[0];
	double h_y = fPeriodicity[1];
	double h_z = fPeriodicity[2];
	
	/* fetch prospective neighbors */
	double target[3];
	target[0] = fCoords(n,0);
	target[1] = fCoords(n,1);
	target[2] = fCoords(n,2);

	/* main cell hits */
	ProcessHitsSorted(target, tol, skiptag, fSortedHits);
	
	/* overlap flags */
	int ix = (h_x < 0.0) ? 0 :
	         ((target[0] + tol > fxmax) ? -1 :
			  ((target[0] - tol < fxmin) ? 1 : 0));
	int iy = (h_y < 0.0) ? 0 :
	         ((target[1] + tol > fymax) ? -1 :
			  ((target[1] - tol < fymin) ? 1 : 0));
	int iz = (h_z < 0.0) ? 0 :
	         ((target[2] + tol > fzmax) ? -1 :
			  ((target[2] - tol < fzmin) ? 1 : 0));

	/* periodic extensions */
	if (ix != 0)
	{
		target[0] += ix*h_x;
		ProcessHitsSorted(target, tol, skiptag, fSortedHits);
		
		if (iy != 0)
		{
			target[1] += iy*h_y;
			ProcessHitsSorted(target, tol, skiptag, fSortedHits);

			if (iz != 0)
			{
				target[2] += iz*h_z;
				ProcessHitsSorted(target, tol, skiptag, fSortedHits);
				target[2] -= iz*h_z;
			}

			target[1] -= iy*h_y;
		}

		if (iz != 0)
		{
			target[2] += iz*h_z;
			ProcessHitsSorted(target, tol, skiptag, fSortedHits);
			target[2] -= iz*h_z;
		}
			
		target[0] -= ix*h_x;
	}

	if (iy != 0)
	{
		target[1] += iy*h_y;
		ProcessHitsSorted(target, tol, skiptag, fSortedHits);
			
		if (iz != 0)
		{
			target[2] += iz*h_z;
			ProcessHitsSorted(target, tol, skiptag, fSortedHits);
			target[2] -= iz*h_z;
		}

		target[1] -= iy*h_y;
	}

	if (iz != 0)
	{
		target[2] += iz*h_z;
		ProcessHitsSorted(target, tol, skiptag, fSortedHits);
		target[2] -= iz*h_z;
	}
	
	/* copy to output */
	fSortedHits.CopyInto(neighbors);
		
	return fSortedHits.Length();
}
