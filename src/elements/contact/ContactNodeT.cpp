/*  $Id: ContactNodeT.cpp,v 1.17 2003-02-03 04:40:18 paklein Exp $ */
#include "ContactNodeT.h"

#include "FaceT.h"

using namespace Tahoe;

ContactNodeT::ContactNodeT(ContactSurfaceT& surface, int node_tag):
	fSurface(surface)
{
	fNodeTag         = node_tag;
	fStatus 	 	= kNoProjection;
	fOpposingSurface = NULL;
	fOpposingFace    = NULL;
	fOriginalOpposingFace    = NULL;
	fxi[0]           = 0.0 ;
	fxi[1]           = 0.0 ;
	fGap             = 1.e8; // NEED TO FIX THIS
	fPressure        = 0.0; 
}

ContactNodeT::~ContactNodeT(void)
{
}

void
ContactNodeT::PrintData(ostream& out)
{
	out << "ContactNode "<< fNodeTag 
	    << " opposing surface " << fOpposingSurface->Tag()
	    << " gap " << fGap 
	    << " xi " << fxi[0] << " " << fxi[1] << '\n';
}

bool
ContactNodeT::AssignOpposing
(const SurfaceT& opposing_surface, const FaceT& opposing_face,
double* xi, double g)
{ // should compare to see if better, (requires initialization)
	fStatus = kProjection;
	/* cast SurfaceT to ContactSurfaceT */
	fOpposingSurface = ((ContactSurfaceT*) &opposing_surface) ;
	fOpposingFace    = &opposing_face ;
	fxi[0] = xi[0] ;
	if (fOpposingSurface->NumSD() == 3 ) {fxi[1] = xi[1] ; }
	fGap = g ;
#if 0
	PrintData(cout);
#endif
	return 1;
}

void 
ContactNodeT::AssignOriginal(void)
{ 
	fOriginalOpposingFace = fOpposingFace;
	fxiO[0] = fxi[0];
	fxiO[1] = fxi[1]; 
}


void 
ContactNodeT::ComputeSlip(double* slip)
{
	
	/* current position of contact point on face */
	if (fOriginalOpposingFace) {
	 double x2_O [3] ;	
	 fOriginalOpposingFace->InterpolatePosition(fxiO,x2_O);
	 /* current position of node */
	 const double* x1 =fSurface.Position(fNodeTag);	
	 slip[0] = x2_O[0] - x1[0];
	 slip[1] = x2_O[1] - x1[1];
	 if (fSurface.NumSD()==3) {slip[2] = x2_O[2] - x1[2];}
	}
	else {
	 slip[0] = 0.0; slip[1] = 0.0;
	 if (fSurface.NumSD()==3) {slip[2] = 0.0;}
	}
}

