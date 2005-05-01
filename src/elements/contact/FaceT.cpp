/* $Id: FaceT.cpp,v 1.12 2005-05-01 20:30:45 paklein Exp $ */
#include "FaceT.h"

#include "SurfaceT.h" // this is for global nodes

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FaceT**>::fByteCopy = true;
} /* namespace Tahoe */

/*constructor*/
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


void
FaceT::InterpolateVector
(const dArrayT& local_coordinates, const dArray2DT& nodal_vectors, 
 dArrayT& vector) const
{
	dArrayT shape_f(nodal_vectors.MajorDim());
	ComputeShapeFunctions (local_coordinates.Pointer(), shape_f);
	vector = 0.0;
	for (int i=0; i<nodal_vectors.MajorDim(); i++) {
		for (int j=0; j<nodal_vectors.MinorDim(); j++) {
			vector[j] += shape_f[i]*nodal_vectors(i,j);
		}
	}
}

