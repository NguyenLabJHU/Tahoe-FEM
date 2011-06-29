/* $Id: ThreeBodyPropertyT_factory.cpp,v 1.1 2004-11-23 01:43:20 cjkimme Exp $ */
#include "ThreeBodyPropertyT.h"
#include <string.h>

/* subclasses supporting the factory method */
#include "StillingerWeberT.h"

using namespace Tahoe;

/* pair property factor */
ThreeBodyPropertyT* ThreeBodyPropertyT::New(const char* name, const BasicSupportT* support)
{
	if (strcmp(name, "Stillinger_Weber") == 0)
		return new StillingerWeberT;
	else
		return NULL;
}
