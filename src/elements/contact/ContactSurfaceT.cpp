/*  $Id: ContactSurfaceT.cpp,v 1.3 2001-04-19 23:47:01 rjones Exp $ */
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
	for(int i = 0; i < fContactNodes.Length(); i++){
		fContactNodes[i] = new ContactNodeT(*this,i);
	}
#if 0
	if (friction)
		fPreviousContactPoints.Allocate(fGlobalNodes.Length());
#endif
}

void
ContactSurfaceT::CopyCurrentToPrevious()
{
	for (int i = 0 ; i < fContactNodes.Length() ; i++) {
#if 0
		fPreviousContactPoints[i] = fContactPoints[i];
		fContactPoints[i].OpposingSurface() = NULL;
#endif
	}
}
