/* $Id: DiffusionMatListT.cpp,v 1.5 2003-06-09 06:53:11 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "DiffusionMatListT.h"
#include "DiffusionMatSupportT.h"
#include "fstreamT.h"

/* diffusion materials */
#include "DiffusionMaterialT.h"
#include "NLDiffusionMaterialT.h"

using namespace Tahoe;

/* constructors */
DiffusionMatListT::	DiffusionMatListT(int length, const DiffusionMatSupportT& support):
	MaterialListT(length),
	fDiffusionMatSupport(support)
{

}

/* read material data from the input stream */
void DiffusionMatListT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "DiffusionMatListT::ReadMaterialData";

	/* read material data */
	for (int i = 0; i < fLength; i++)
	{
		int matnum, matcode;
		in >> matnum; matnum--;
		in >> matcode;
		
		/* checks */
		if (matnum < 0  || matnum >= fLength) ExceptionT::BadInputValue(caller);

		/* repeated material number */
		if (fArray[matnum] != NULL)
			ExceptionT::BadInputValue(caller, "repeated material number: %d", matnum+1);
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kLinear:
			{
				fArray[matnum] = new DiffusionMaterialT(in, fDiffusionMatSupport);
				break;
			}
			case kNonLinear:
			{
				fArray[matnum] = new NLDiffusionMaterialT(in, fDiffusionMatSupport);
				break;
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
		}

		/* verify construction */
		if (!fArray[matnum]) ExceptionT::OutOfMemory(caller);
	}
}
