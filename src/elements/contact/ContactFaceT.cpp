/* $Id $ */

#include "ContactFaceT.h"

/*constructor*/

using namespace Tahoe;

ContactFaceT::ContactFaceT
(FaceT* face):
	fFace(face)
{
    fMultiplierConnectivity.Dimension(fFace->NumNodes());
}

/*destructor*/ 
ContactFaceT::~ContactFaceT (void) { }
