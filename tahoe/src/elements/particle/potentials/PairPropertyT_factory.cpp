/* $Id: PairPropertyT_factory.cpp,v 1.1.2.1 2004-06-14 04:56:38 paklein Exp $ */
#include "PairPropertyT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "HarmonicPairT.h"
#include "LennardJonesPairT.h"
#include "ParadynPairT.h"
#include "MatsuiPairT.h"

using namespace Tahoe;

/* pair property factor */
PairPropertyT* PairPropertyT::New(const char* name, const BasicSupportT* support)
{
	if (strcmp(name, "harmonic") == 0)
		return new HarmonicPairT;
	else if (strcmp(name, "Lennard_Jones") == 0)
		return new LennardJonesPairT;
	else if (strcmp(name, "Paradyn_pair") == 0)
	{
		/* requires support */
		if (!support) ExceptionT::GeneralFail("PairPropertyT::New", "\"%s\" requires BasicSupportT", name);
		
		return new ParadynPairT(*support);	
	}
	else if (strcmp(name, "Matsui") == 0)
		return new MatsuiPairT;
	else
		return NULL;
}
