/* $Id: ThreeBodyPropertyT_factory.cpp,v 1.2 2011-12-01 21:11:39 bcyansfn Exp $ */
#include "ThreeBodyPropertyT.h"
#include <cstring>

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
