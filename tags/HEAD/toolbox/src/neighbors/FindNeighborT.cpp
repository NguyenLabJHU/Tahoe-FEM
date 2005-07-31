/* $Id: FindNeighborT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (03/21/1997)                                          */
/* FindNeighborT.cpp                                                      */

#include "FindNeighborT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"
#include "iGridManager2DT.h"
#include "iGridManager3DT.h"

/* search grid parameters */
const int Grid_x_count = 20;
const int Grid_y_count = 20;
const int Grid_z_count = 20;

/* Constructor */
FindNeighborT::FindNeighborT(const dArray2DT& coords, int maxneighbors):
	fMaxNeighbors(maxneighbors),
	fNodeMap(NULL),
	fglCoords(coords),
	fGrid2D(NULL),
	fGrid3D(NULL)
{
	if (fMaxNeighbors < 1) throw eBadInputValue;

	fCoords.Alias(fglCoords);
	Allocate(fCoords.MajorDim(), fCoords.MinorDim());
			
	/* set search grid */
	if (fnsd == 2)
	{
		fGrid2D = new iGridManager2DT(Grid_x_count, Grid_y_count, fCoords, NULL);	
		if (!fGrid2D) throw eOutOfMemory;		
	}
	else
	{
		fGrid3D = new iGridManager3DT(Grid_x_count, Grid_y_count, Grid_z_count, fCoords, NULL);
		if (!fGrid3D) throw eOutOfMemory;		
	}
}

FindNeighborT::FindNeighborT(const iArrayT& nodesused, const dArray2DT& coords,
	int maxneighbors):
	fMaxNeighbors(maxneighbors),
	fNodeMap(&nodesused),
	fglCoords(coords),
	fGrid2D(NULL),
	fGrid3D(NULL)
{
	if (fMaxNeighbors < 1) throw eBadInputValue;

	Allocate(fNodeMap->Length(), coords.MinorDim());
	
	/* collect coords of used nodes */
	fCoords.Allocate(fNodeMap->Length(), fglCoords.MinorDim());
	fCoords.RowCollect(*fNodeMap, fglCoords);

	/* set search grid */
	if (fnsd == 2)
	{
		fGrid2D = new iGridManager2DT(Grid_x_count, Grid_y_count, fCoords, 0);	
		if (!fGrid2D) throw eOutOfMemory;		
	}
	else
	{
		fGrid3D = new iGridManager3DT(Grid_x_count, Grid_y_count, Grid_z_count, fCoords, 0);	
		if (!fGrid3D) throw eOutOfMemory;		
	}
}

/* Destructor */
FindNeighborT::~FindNeighborT(void)
{
	delete fGrid2D;
	delete fGrid3D;
}

/* Print neighbors to output stream */
void FindNeighborT::OutputNeighors(ostream& out, double tolerance)
{
	/* find neighbors */
	if (fnsd == 2)
		FindNeighors2D(tolerance);
	else
		FindNeighors3D(tolerance);
	
	/* will not print duplicates */
	int count = 0;
	for (int node = 0; node < fNumPts; node++)
		for (int neigh = 0; neigh < fCount[node]; neigh++)
			if (fNeighbors(node, neigh) > node)
			{
				out << setw(kIntWidth) << ++count;
				out << setw(kIntWidth) << node+1;
				out << setw(kIntWidth) << fNeighbors(node, neigh) << '\n';
			}
			
	out << '\n';
}

void FindNeighborT::GetNeighors(iArray2DT& edges, double tolerance)
{
	/* find neighbors - exhaustive search */
	if (fnsd == 2)
		FindNeighors2D(tolerance);
	else
		FindNeighors3D(tolerance);
	
	/* copy unique edges */
	if (fNodeMap)
	{
		edges.Allocate(Count2BodyMapped(), 2);
		Set2BodyMapped(edges);
	}
	else
	{
		edges.Allocate(Count2Body(), 2);	
		Set2Body(edges);
	}
}

/* Print coordinate data to the output stream */
void  FindNeighborT::PrintCoords(ostream& out) const
{
	int d_width = OutputWidth(out, fCoords.Pointer());
	for (int node = 0; node < fNumPts; node++)
	{
		out << setw(kIntWidth) << node+1;

		for (int nsd = 0; nsd < fnsd; nsd++)
			out << setw(d_width) << fCoords(nsd, node);
			
		out << '\n';
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* Determine neighbors */
void FindNeighborT::FindNeighors2D(double tolerance)
{
	/* update search coords */
	if (fNodeMap) fCoords.RowCollect(*fNodeMap, fglCoords);

	/* reset the grid */
	fGrid2D->Reset();

	/* using the grid */
	AutoArrayT<int> neighbors;
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		fGrid2D->Neighbors(i, tolerance, neighbors);
		fCount[i] = neighbors.Length();
		
		if (neighbors.Length() > fNeighbors.MinorDim()) throw eSizeMismatch;
		
		/* copy in */
		memcpy(fNeighbors(i), neighbors.Pointer(), sizeof(int)*neighbors.Length());
	}

	/* check that all are connected */
	int noneighbors = fCount.Count(0);
	if (noneighbors > 0)
	{
		cout << "\nFindNeighborT::FindNeighors2D: " << noneighbors;
		cout << " without neighnors" << endl;
	}
}

void FindNeighborT::FindNeighors3D(double tolerance)
{
	/* update search coords */
	if (fNodeMap) fCoords.RowCollect(*fNodeMap, fglCoords);

	/* reset the grid */
	fGrid3D->Reset();

	/* using the grid */
	AutoArrayT<int> neighbors;
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		fGrid3D->Neighbors(i, tolerance, neighbors);
		fCount[i] = neighbors.Length();
		
		if (neighbors.Length() > fNeighbors.MinorDim()) throw eSizeMismatch;
		
		/* copy in */
		memcpy(fNeighbors(i), neighbors.Pointer(), sizeof(int)*neighbors.Length());
	}

	/* check that all are connected */
	int noneighbors = fCount.Count(0);
	if (noneighbors > 0)
	{
		cout << "\nFindNeighborT::FindNeighors3D: " << noneighbors;
		cout << " without neighnors" << endl;
	}
}

/* Allocate memory */
void FindNeighborT::Allocate(int numpts, int nsd)
{
	fNumPts = numpts;
	fnsd    = nsd;
	
	/* only support 2D and 3D data */
	if (fnsd != 2 && fnsd != 3) throw eBadInputValue;

	/* Allocate memory */
	fCount.Allocate(fNumPts);
	fNeighbors.Allocate(fNumPts, fMaxNeighbors);
}

/* Determine number of 2 body interactions */
int FindNeighborT::Count2Body(void) const
{
	int count = 0;
	for (int node = 0; node < fNumPts; node++)
	{
		int* pneigh = fNeighbors(node);	
		int  num    = node;
	
		for (int neigh = 0; neigh < fCount[node]; neigh++)
			if (*pneigh++ > num) count++;
	}
				
	return count;
}

/* Determine number of 2 body interactions */
int FindNeighborT::Count2BodyMapped(void) const
{
	int count = 0;
	for (int node = 0; node < fNumPts; node++)
	{
		int* pneigh = fNeighbors(node);
		int  num    = (*fNodeMap)[node];
	
		for (int neigh = 0; neigh < fCount[node]; neigh++)
			if ((*fNodeMap)[*pneigh++] > num) count++;
	}
				
	return count;
}

/* copy unique 2 body into edges */
void FindNeighborT::Set2Body(iArray2DT& edges) const
{
	int* pedges = edges.Pointer();

	for (int node = 0; node < fNumPts; node++)
	{
		int* pneigh = fNeighbors(node);
	
		for (int neigh = 0; neigh < fCount[node]; neigh++)
		{		
			if (*pneigh > node)
			{
				pedges[0] = node;
				pedges[1] = *pneigh;
			
				pedges += 2;
			}

			pneigh++;
		}
	}
}

void FindNeighborT::Set2BodyMapped(iArray2DT& edges) const
{
	int* pedges = edges.Pointer();

	for (int node = 0; node < fNumPts; node++)
	{
		int* pneigh = fNeighbors(node);
		
		for (int neigh = 0; neigh < fCount[node]; neigh++)
		{
			if ((*fNodeMap)[*pneigh] > (*fNodeMap)[node])
			{
				pedges[0] = (*fNodeMap)[node];
				pedges[1] = (*fNodeMap)[*pneigh];
			
				pedges += 2;
			}
			
			pneigh++;
		}
	}
}
