/* $Id: LocalArrayGroupT.cpp,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (09/11/1998)                                          */
/* Class to manage a list of group of dynamically changing size           */
/* LocalArrayT's                                                          */
/* NOTE: all registered arrays will be shallow.                           */

#include "LocalArrayGroupT.h"
#include "LocalArrayT.h"

/* constructor */
LocalArrayGroupT::LocalArrayGroupT(int headroom):
	MemoryGroupT<double>(headroom),
	fNumNodes(0),
	fMinorDim(0)
{

}

/* add array to list of managed */
void LocalArrayGroupT::Register(LocalArrayT& localarray)
{
	/* must all have the same minor dimension */
	if (fArrays.Length() > 0)
	{
		/* size check */
		if (fMinorDim != localarray.MinorDim())
		{
			cout << "\n LocalArrayGroupT::Register: all arrays must";
			cout << " be of the same minor dimension: ";
			cout << fMinorDim << endl;
			throw eSizeMismatch;
		}
	}
	else
		fMinorDim = localarray.MinorDim();

	/* inherited */
	MemoryGroupT<double>::Register(localarray);
}

/* set number of nodes */
void LocalArrayGroupT::SetNumberOfNodes(int numnodes)
{
	if (numnodes != fNumNodes)
	{
		fNumNodes = numnodes;
	
		/* need more memory */
		int blocksize = fNumNodes*fMinorDim;
		if (blocksize > BlockSize()) SetBlockSize(blocksize, false);
		
		/* reset dimensions */
		for (int i = 0; i < fArrays.Length(); i++)
		{
			/* safe cast due to type filtering by Register */
			LocalArrayT* parray = (LocalArrayT*) fArrays[i];
		
			/* set internal parameters */
			parray->Set(fNumNodes, fMinorDim, BlockPointer(i));
		}
	}
}
