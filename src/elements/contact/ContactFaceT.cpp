/* $Id $ */

#include "ContactFaceT.h"

/*constructor*/
ContactFaceT::ContactFaceT
(FaceT* face):
	fFace(face)
{
    fMultiplierConnectivity.Allocate(fFace->NumNodes());
}

/*destructor*/ 
ContactFaceT::~ContactFaceT (void) { }
