/* $Id: FaceT.cpp,v 1.7 2002-10-20 22:48:21 paklein Exp $ */

#include "FaceT.h"

#include "SurfaceT.h" // this is for global nodes

/*constructor*/

using namespace Tahoe;

FaceT::FaceT
(SurfaceT& surface,dArray2DT& surface_coordinates,	
int num_face_nodes, int* connectivity):
	fSurface(surface),
	fSurfaceCoordinates(surface_coordinates)
{
        fConnectivity.Dimension(num_face_nodes);
        fGlobalConnectivity.Dimension(num_face_nodes);
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
