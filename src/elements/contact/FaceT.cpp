/* $Id: FaceT.cpp,v 1.2 2001-04-09 22:28:55 rjones Exp $ */

#include "FaceT.h"

/*constructor*/
FaceT::FaceT
(SurfaceT& surface,dArray2DT& surface_coordinates,	
int num_face_nodes, int* connectivity):
	fSurface(surface),
	fSurfaceCoordinates(surface_coordinates)
{
        fNumNodes = num_face_nodes;
        fConnectivity.Allocate(fNumNodes);
        for (int i = 0; i < fNumNodes; i++) {
                fConnectivity[i] = connectivity[i];
	}
}

/*destructor*/ 
FaceT::~FaceT (void) { }
