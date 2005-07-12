/* $Id: ThreeBodyPropertyT_factory.cpp,v 1.1.14.1 2005-06-10 23:02:32 paklein Exp $ */
#include "ThreeBodyPropertyT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "StillingerWeberT.h"

using namespace Tahoe;

/* pair property factor */
ThreeBodyPropertyT* ThreeBodyPropertyT::New(const char* name, const BasicSupportT* support)
{
#pragma unused(support)
	if (strcmp(name, "Stillinger_Weber") == 0)
		return new StillingerWeberT;
	else
		return NULL;
}
