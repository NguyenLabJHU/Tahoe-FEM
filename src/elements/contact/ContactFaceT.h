/* $Id: ContactFaceT.h,v 1.2 2002-07-02 19:55:19 cjkimme Exp $ */

#ifndef _CONTACT_FACE_T_H_
#define _CONTACT_FACE_T_H_

/* direct members */
#include "iArrayT.h"
#include "ContactSurfaceT.h"

/* forward declarations */

namespace Tahoe {

class FaceT;
class dMatrixT;

class ContactFaceT 
{
public:

	/* constructor */
	ContactFaceT (FaceT* face); 	

	/* (virtual) destructor */
  	virtual ~ContactFaceT(void);

	virtual void ComputePressureFunctions
		(const double* local_coordinates, dMatrixT& shape_functions)
        const=0;

	inline const iArrayT& MultiplierConnectivity(void) const 
		{return fMultiplierConnectivity;}
	inline void SetMultiplierConnectivity(void) 
		{ContactSurfaceT& surface = (ContactSurfaceT&) fFace->Surface();
		surface.MultiplierTags(fFace->Connectivity(),fMultiplierConnectivity);}

protected:
	FaceT* fFace;
	iArrayT fMultiplierConnectivity;
};

} // namespace Tahoe 
#endif /* _CONTACT_FACE_T_H_ */

