/*  $Id: ContactNodeT.cpp,v 1.4 2001-04-24 00:33:21 rjones Exp $ */
#include "ContactNodeT.h"

#include "SurfaceT.h"
#include "FaceT.h"

/* parameters */

ContactNodeT::ContactNodeT(SurfaceT& surface, int node_tag):
	fSurface(surface)
{
	fNodeTag         = node_tag;
	fOpposingSurface = NULL;
	fOpposingFace    = NULL;
	fxi[0]           = 0.0 ;
	fxi[1]           = 0.0 ;
	fGap             = 1.0e8 ; // this needs to be TOL_G?
}

ContactNodeT::~ContactNodeT(void)
{
}

void
ContactNodeT::PrintData(ostream& out)
{
	out << "gap " << fGap << '\n';
	out << "xi " << fxi[0] << " " << fxi[1] << '\n';
}

bool
ContactNodeT::AssignOpposing
(const SurfaceT* opposing_surface, const FaceT* opposing_face,
double* xi, double g)
{
	if ( g < fGap ) { // this is a little naive
        	fOpposingSurface = opposing_surface ;
        	fOpposingFace    = opposing_face ;
        	fxi[0] = xi[0] ;
		if (fOpposingSurface->NumSD() == 3 ) 
			{fxi[1] = xi[1] ; }
        	fGap = g ;
		PrintData(cout);
		return 1;
	}
	return 0;
}

void 
ContactNodeT::UpdateOpposing
(double* xi, double g)
{
                fxi[0] = xi[0] ;
                if (fOpposingSurface->NumSD() == 3 )
                        {fxi[1] = xi[1] ; }
                fGap = g ;
}

