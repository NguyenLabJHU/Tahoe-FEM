/* $Id $ */

#include "ContactFaceT.h"

/*constructor*/

using namespace Tahoe;

ContactFaceT::ContactFaceT
(FaceT* face):
	fFace(face)
{
    fMultiplierConnectivity.Allocate(fFace->NumNodes());
}

/*destructor*/ 
ContactFaceT::~ContactFaceT (void) { }
