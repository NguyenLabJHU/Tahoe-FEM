/* $Id: LocalArrayGroupT.cpp,v 1.4 2003-01-27 06:42:44 paklein Exp $ */
/* created: paklein (09/11/1998) */
#include "LocalArrayGroupT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* constructor */
LocalArrayGroupT::LocalArrayGroupT(int headroom):
	MemoryGroupT<double>(headroom, true), /* use pooled memory */
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
			ExceptionT::SizeMismatch("LocalArrayGroupT::Register", 
				"all arrays must be of the same minor dimension: %d", fMinorDim);
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
