/* $Id: FindNeighbor23T.cpp,v 1.4 2003-11-21 22:41:59 paklein Exp $ */
/* created: paklein (03/21/1997)                                          */
/* FindNeighbor23T.cpp                                                    */

#include "FindNeighbor23T.h"
#include <iostream.h>

/* Constructors */

using namespace Tahoe;

FindNeighbor23T::FindNeighbor23T(const dArray2DT& coords, int maxneighbors):
	FindNeighborT(coords, maxneighbors)
{

}

FindNeighbor23T::FindNeighbor23T(const iArrayT& nodesused,
	const dArray2DT& coords, int maxneighbors):
	FindNeighborT(nodesused, coords, maxneighbors)
{

}

/* Print neighbors to output stream */
void FindNeighbor23T::OutputNeighors(ostream& out, double tolerance)
{
	/* find and print 2 body */
	FindNeighborT::OutputNeighors(out, tolerance);

	/* count number of 3 body and allocate */
	iArray2DT angles(Count3Body(), 3);
	
	/* determine 3 body lists */
	Set3Body(angles);
	
	/* output */
	angles.WriteNumbered(out);
}

void FindNeighbor23T::GetNeighors(iArray2DT& edges, iArray2DT& angles,
	double tolerance)
{
	/* find and print 2 body */
	FindNeighborT::GetNeighors(edges, tolerance);

	/* count number of 3 body and allocate */
	angles.Dimension(Count3Body(), 3);
	
	/* determine 3 body lists */
	if (fNodeMap)
		Set3BodyMapped(angles);
	else
		Set3Body(angles);
}	

/**********************************************************************
* Private
**********************************************************************/

/* Determine number of 3 body interactions */
int FindNeighbor23T::Count3Body(void) const
{
	/* defines f(n) = 1 + 2 + ... + n-1 */
	iArrayT f(fMaxNeighbors + 1); //OFFSET?
	
	f[0] = 0;
	for (int i = 1; i <= fMaxNeighbors; i++)
		f[i] = f[i-1] + i-1;

	int count = 0;	
	for (int j = 0; j < fNumPts; j++)
		count += f[ fCount[j] ];

	return count;
}

/* Determine 3 body interactions */
void FindNeighbor23T::Set3Body(iArray2DT& angles) const
{
	int* p3Body = angles(0);
	for (int node = 0; node < fNumPts; node++)
	{
		const int* pneigh = fNeighbors(node);
	
		for (int i = 0; i < fCount[node]; i++)
			for (int j = i + 1; j < fCount[node]; j++)
			{	
				p3Body[0] = pneigh[i];
				p3Body[1] = node;
				p3Body[2] = pneigh[j];
				
				p3Body += 3;
			}
	}
}

/* Determine 3 body interactions */
void FindNeighbor23T::Set3BodyMapped(iArray2DT& angles) const
{
	int* p3Body = angles(0);
	for (int node = 0; node < fNumPts; node++)
	{
		const int* pneigh = fNeighbors(node);
	
		for (int i = 0; i < fCount[node]; i++)
			for (int j = i + 1; j < fCount[node]; j++)
			{	
				p3Body[0] = (*fNodeMap)[pneigh[i]];
				p3Body[1] = (*fNodeMap)[node];
				p3Body[2] = (*fNodeMap)[pneigh[j]];
				
				p3Body += 3;
			}
	}
}
