/* $Id: iPeriodicGrid3DT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/18/1997)                                          */

#include "iPeriodicGrid3DT.h"
#include "dArrayT.h"
#include "iArrayT.h"

/* constructor */
iPeriodicGrid3DT::iPeriodicGrid3DT(int nx, int ny, int nz,
	const dArray2DT& coords, const iArrayT* nodes_used,
	const dArrayT& periodicity):
	iGridManager3DT(nx, ny, nz, coords, nodes_used),
	fPeriodicity(periodicity)
{

}	

/* neighbors - returns the number of neighbors */
int iPeriodicGrid3DT::PeriodicNeighbors(int n, double tol, iArrayT& neighbors)
{
	/* initialize */
	fSortedHits.Allocate(0);
	int skiptag = n;

	double hx = fPeriodicity[0];
	double hy = fPeriodicity[1];
	double hz = fPeriodicity[2];
	
	/* fetch prospective neighbors */
	double* target = fCoords(n);

	/* main cell hits */
	ProcessHitsSorted(target, tol, skiptag, fSortedHits);
	
	/* overlap flags */
	int ix = (hx < 0.0) ? 0 :
	         ((target[0] + tol > fxmax) ? -1 :
			  ((target[0] - tol < fxmin) ? 1 : 0));
	int iy = (hy < 0.0) ? 0 :
	         ((target[1] + tol > fymax) ? -1 :
			  ((target[1] - tol < fymin) ? 1 : 0));
	int iz = (hz < 0.0) ? 0 :
	         ((target[2] + tol > fzmax) ? -1 :
			  ((target[2] - tol < fzmin) ? 1 : 0));

	/* periodic extensions */
	if (ix != 0)
	{
		target[0] += ix*hx;
		ProcessHitsSorted(target, tol, skiptag, fSortedHits);
		
		if (iy != 0)
		{
			target[1] += iy*hy;
			ProcessHitsSorted(target, tol, skiptag, fSortedHits);

			if (iz != 0)
			{
				target[2] += iz*hz;
				ProcessHitsSorted(target, tol, skiptag, fSortedHits);
				target[2] -= iz*hz;
			}

			target[1] -= iy*hy;
		}

		if (iz != 0)
		{
			target[2] += iz*hz;
			ProcessHitsSorted(target, tol, skiptag, fSortedHits);
			target[2] -= iz*hz;
		}
			
		target[0] -= ix*hx;
	}

	if (iy != 0)
	{
		target[1] += iy*hy;
		ProcessHitsSorted(target, tol, skiptag, fSortedHits);
			
		if (iz != 0)
		{
			target[2] += iz*hz;
			ProcessHitsSorted(target, tol, skiptag, fSortedHits);
			target[2] -= iz*hz;
		}

		target[1] -= iy*hy;
	}

	if (iz != 0)
	{
		target[2] += iz*hz;
		ProcessHitsSorted(target, tol, skiptag, fSortedHits);
		target[2] -= iz*hz;
	}
	
	/* copy to output */
	fSortedHits.CopyInto(neighbors);
		
	return fSortedHits.Length();
}
