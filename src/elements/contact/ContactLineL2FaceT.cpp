/* $Id: ContactLineL2FaceT.cpp,v 1.2 2002-07-02 19:55:19 cjkimme Exp $ */

#include "ContactLineL2FaceT.h"

#include "dMatrixT.h"


using namespace Tahoe;

ContactLineL2FaceT::ContactLineL2FaceT
(FaceT* face):
	ContactFaceT(face)
{
}

ContactLineL2FaceT::~ContactLineL2FaceT (void)
{
}

void
ContactLineL2FaceT::ComputePressureFunctions
(const double* local_coordinates, dMatrixT& shape_functions) const
{
	shape_functions = 0.0;
	double xi  = local_coordinates[0];
	shape_functions(0,0) = 0.5 * (1.0 - xi );
	shape_functions(1,0) = 0.5 * (1.0 + xi );
}

