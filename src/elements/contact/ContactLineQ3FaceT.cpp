/* $Id: ContactLineQ3FaceT.cpp,v 1.1 2001-09-24 20:43:25 rjones Exp $ */

#include "ContactLineQ3FaceT.h"

#include "dMatrixT.h"

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

