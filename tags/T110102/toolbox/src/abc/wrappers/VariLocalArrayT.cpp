/* $Id: VariLocalArrayT.cpp,v 1.3 2002-10-20 22:38:57 paklein Exp $ */
/* created: paklein (04/26/1999)                                          */
/* Wrapper for dynamically re-sizing the number of nodes in               */
/* a LocalArrayT's.                                                       */

#include "VariLocalArrayT.h"
#include "LocalArrayT.h"

/* constructors */

using namespace Tahoe;

VariLocalArrayT::VariLocalArrayT(void):
	fMinorDim(0),
	fWard(NULL)
{

}

VariLocalArrayT::VariLocalArrayT(int headroom, LocalArrayT& ward,
	int minordim): fWard(NULL)
{
	SetWard(headroom, ward, minordim);
}

void VariLocalArrayT::SetWard(int headroom, LocalArrayT& ward,
	int minordim)
{
	/* inherited */
	SetHeadRoom(headroom);

	/* can only be called once */
	if (!fWard)
	{
		fMinorDim = minordim;
		fWard     = &ward;
		if (fWard->MinorDim() > 0)
		{
			/* consistency check */
			if (fWard->MinorDim() != fMinorDim) throw ExceptionT::kSizeMismatch;
		}
		else
			/* set minor dimension */
			fWard->Set(0, fMinorDim, NULL);
	}
	else
		throw ExceptionT::kGeneralFail;
}

/* set number of nodes */
void VariLocalArrayT::SetNumberOfNodes(int numnodes)
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	if (numnodes != fWard->NumberOfNodes())
	{
		/* update ArrayT data */
		SetAlias(*fWard, numnodes*fMinorDim, false);
	
		/* update LocalArrayT data */
		fWard->Set(numnodes, fMinorDim, fWard->Pointer());
	}
}
