/* $Id: FaceT.cpp,v 1.4 2001-04-19 23:47:01 rjones Exp $ */

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
	fnormal[0] = 0.0;
	fnormal[1] = 0.0;
	fnormal[2] = 0.0;
}

/*destructor*/ 
FaceT::~FaceT (void) { }
