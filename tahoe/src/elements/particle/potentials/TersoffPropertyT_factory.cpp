/* $Id: TersoffPropertyT_factory.cpp,v 1.1.2.1 2006-06-12 18:40:06 d-farrell2 Exp $ */
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
