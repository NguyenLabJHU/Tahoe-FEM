/* $Id: TersoffPropertyT_factory.cpp,v 1.2.4.1 2011-10-30 06:26:11 bcyansfn Exp $ */
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
