/* $Id: PairPropertyT_factory.cpp,v 1.3 2011-12-01 21:11:39 bcyansfn Exp $ */
#include "PairPropertyT.h"
#include <cstring>

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
		return new ParadynPairT(support);	
	else if (strcmp(name, "Matsui") == 0)
		return new MatsuiPairT;
	else
		return NULL;
}
