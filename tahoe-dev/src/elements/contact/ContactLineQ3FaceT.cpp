/* $Id: ContactLineQ3FaceT.cpp,v 1.2 2002-07-02 19:55:19 cjkimme Exp $ */

#include "ContactLineQ3FaceT.h"

#include "dMatrixT.h"


using namespace Tahoe;

ContactLineQ3FaceT::ContactLineQ3FaceT
(FaceT* face):
	ContactFaceT(face)
{
}

ContactLineQ3FaceT::~ContactLineQ3FaceT (void)
{
}

void
ContactLineQ3FaceT::ComputePressureFunctions
(const double* local_coordinates, dMatrixT& shape_functions) const
{
	shape_functions = 0.0;
    double xi  = local_coordinates[0];
    shape_functions(0,0) = 0.5 * xi * ( xi - 1.0 );
    shape_functions(1,0) = 0.5 * xi * ( xi + 1.0 );
    shape_functions(2,0) = 1.0 - xi * xi ;
}

