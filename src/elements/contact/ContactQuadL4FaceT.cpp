/* $Id: ContactQuadL4FaceT.cpp,v 1.1 2001-09-24 20:43:25 rjones Exp $ */

#include "ContactQuadL4FaceT.h"

#include "dMatrixT.h"

ContactQuadL4FaceT::ContactQuadL4FaceT
(FaceT* face):
	ContactFaceT(face)
{
}

ContactQuadL4FaceT::~ContactQuadL4FaceT (void)
{
}

void
ContactQuadL4FaceT::ComputePressureFunctions
(const double* local_coordinates, dMatrixT& shape_functions) const
{
	shape_functions = 0.0;
    double xi  = local_coordinates[0];
    double eta = local_coordinates[1];
    shape_functions(0,0) = 0.25 * (1.0 - xi ) * (1.0 - eta) ;
    shape_functions(1,0) = 0.25 * (1.0 + xi ) * (1.0 - eta) ;
    shape_functions(2,0) = 0.25 * (1.0 + xi ) * (1.0 + eta) ;
    shape_functions(3,0) = 0.25 * (1.0 - xi ) * (1.0 + eta) ;

}

