/* $Id: FaceT.cpp,v 1.5 2001-05-23 14:45:04 rjones Exp $ */

#include "FaceT.h"

#include "SurfaceT.h" // this is for global nodes

/*constructor*/
FaceT::FaceT
(SurfaceT& surface,dArray2DT& surface_coordinates,	
int num_face_nodes, int* connectivity):
	fSurface(surface),
	fSurfaceCoordinates(surface_coordinates)
{
        fConnectivity.Allocate(num_face_nodes);
        fGlobalConnectivity.Allocate(num_face_nodes);
        for (int i = 0; i < num_face_nodes; i++) {
                fConnectivity[i] = connectivity[i];
                // NOTE : this is ugly
                fGlobalConnectivity[i] 
		   = surface.GlobalNodes()[connectivity[i]];
	}
	fnormal[0] = 0.0;
	fnormal[1] = 0.0;
	fnormal[2] = 0.0;
}

/*destructor*/ 
FaceT::~FaceT (void) { }
