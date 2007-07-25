/* $Id: TersoffPropertyT_factory.cpp,v 1.2 2006-07-25 16:29:47 d-farrell2 Exp $ */
#include "TersoffPropertyT.h"
#include <string.h>

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
