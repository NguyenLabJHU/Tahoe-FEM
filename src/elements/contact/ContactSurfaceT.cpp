/*  $Id: ContactSurfaceT.cpp,v 1.2 2001-04-16 17:30:51 rjones Exp $ */
#include "ContactSurfaceT.h"

#include "SurfaceT.h"
#include "ContactNodeT.h"

/* parameters */

ContactSurfaceT::ContactSurfaceT(void)
{
}

ContactSurfaceT::~ContactSurfaceT(void)
{
}

void
ContactSurfaceT::AllocateContactNodes()
{
	fContactNodes.Allocate(fGlobalNodes.Length());
}
