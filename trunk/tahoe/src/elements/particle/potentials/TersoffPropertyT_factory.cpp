/* $Id: TersoffPropertyT_factory.cpp,v 1.3 2011-12-01 21:11:39 bcyansfn Exp $ */
#include "TersoffPropertyT.h"
#include <cstring>

/* subclasses supporting the factory method */
#include "TersoffPairT.h"

using namespace Tahoe;

/* pair property factor */
TersoffPropertyT* TersoffPropertyT::New(const char* name, const BasicSupportT* support)
{
	if (strcmp(name, "Tersoff") == 0)
		return new TersoffPairT;
	else
		return NULL;
}
