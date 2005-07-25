/* $Id: ThreeBodyPropertyT_factory.cpp,v 1.1.24.1 2005-07-25 02:37:14 paklein Exp $ */
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
