/* $Id: TiedPotentialBaseT.cpp,v 1.1 2003-03-26 20:00:59 cjkimme Exp $  */
/* created: cjkimme (10/23/2001) */

#include "TiedPotentialBaseT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "iArrayT.h"

using namespace Tahoe;

/* constructor */
TiedPotentialBaseT::TiedPotentialBaseT(void):
	iBulkGroups(0)
{
	// Do nothing
}

/* destructor */
TiedPotentialBaseT::~TiedPotentialBaseT(void) 
{
	// Nope 
}


iArrayT& TiedPotentialBaseT::BulkGroups(void)
{
  	return iBulkGroups;
}


