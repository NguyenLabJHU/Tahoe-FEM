/* $Id: FaceT.cpp,v 1.3 2001-04-11 14:48:57 rjones Exp $ */

#include "FaceT.h"

/*constructor*/
FaceT::FaceT
(SurfaceT& surface,dArray2DT& surface_coordinates,	
int num_face_nodes, int* connectivity):
	fSurface(surface),
	fSurfaceCoordinates(surface_coordinates)
{
        fConnectivity.Allocate(num_face_nodes);
        for (int i = 0; i < num_face_nodes; i++) {
                fConnectivity[i] = connectivity[i];
	}
}

/*destructor*/ 
FaceT::~FaceT (void) { }
