/*  $Id: ContactNodeT.cpp,v 1.12 2001-09-10 23:26:18 rjones Exp $ */
#include "ContactNodeT.h"

#include "FaceT.h"
#include "ContactElementT.h"

/* parameters */

ContactNodeT::ContactNodeT(SurfaceT& surface, int node_tag):
	fSurface(surface)
{
	fNodeTag         = node_tag;
	fStatus 	 = kNoProjection;
	fOpposingSurface = NULL;
	fOpposingFace    = NULL;
	fOriginalOpposingFace    = NULL;
	fxi[0]           = 0.0 ;
	fxi[1]           = 0.0 ;
	fGap             = 1.e8; // NEED TO FIX THIS
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
	fStatus = kContact;
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
ContactNodeT::UpdateOpposing
(double* xi, double g)
{
                fxi[0] = xi[0] ;
                if (fOpposingSurface->NumSD() == 3 )
                        {fxi[1] = xi[1] ; }
                fGap = g ;
}

void 
ContactNodeT::AssignOriginal(void)
{ 
	fOriginalOpposingFace = fOpposingFace;
	fxiO[0] = fxi[0];
	fxiO[1] = fxi[1]; 
}


void 
ContactNodeT::AssignStatus(nMatrixT<dArrayT>& enforcement_parameters)
{
	if (fOpposingSurface) {
		dArrayT& parameters = enforcement_parameters
			(fSurface.Tag(),fOpposingSurface->Tag()) ;
#if 0
		if(fGap < parameters[ContactElementT::ktol_gap]) {
			fStatus = kContact;
		}
		else {
			fStatus = kProjection;
		}
#endif 
		fStatus = kProjection;
	}
	else {
		fStatus = kNoProjection;
	}
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

