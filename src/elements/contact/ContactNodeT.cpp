/*  $Id: ContactNodeT.cpp,v 1.1 2001-04-16 17:30:50 rjones Exp $ */
#include "ContactNodeT.h"

#include "ContactSurfaceT.h"
#include "FaceT.h"

/* parameters */

ContactNodeT::ContactNodeT(void)
{
	fOpposingSurface = NULL;
	fOpposingFace    = NULL;
	fxi[0]           = 0.0 ;
	fxi[1]           = 0.0 ;
	fGap             = 0.0 ; // this needs to be TOL_G
}

ContactNodeT::~ContactNodeT(void)
{
}

void
ContactNodeT::PrintData(ostream& out)
{
	out << "gap " << fGap << '\n';
}

bool
ContactNodeT::AssignOpposing
(SurfaceT* opposing_surface, FaceT* opposing_face,double* xi, double g)
{
	if ( g < fGap ) { // this is a little naive
        	fOpposingSurface = opposing_surface ;
        	fOpposingFace    = opposing_face ;
        	fxi[0] = xi[0] ;
		if (fOpposingSurface->NumSD() == 3 ) 
			{fxi[1] = xi[1] ; }
        	fGap = g ;
		return 1;
	}
	return 0;
}
