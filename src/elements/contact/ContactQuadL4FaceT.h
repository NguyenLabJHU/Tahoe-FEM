/* $Id: ContactQuadL4FaceT.h,v 1.2 2002-07-02 19:55:19 cjkimme Exp $ */

#ifndef _CONTACT_QUADL4_FACE_T_H_
#define _CONTACT_QUADL4_FACE_T_H_

/* base class */
#include "ContactFaceT.h"


namespace Tahoe {

class ContactQuadL4FaceT : public ContactFaceT
{
public:

	/* constructor */
	ContactQuadL4FaceT (FaceT* face); 	

	/* destructor */
  	~ContactQuadL4FaceT(void);

	void ComputePressureFunctions
		(const double* local_coordinates, dMatrixT& shape_functions)
        const;

};

} // namespace Tahoe 
#endif /* _CONTACT_QUADL4_FACE_T_H_ */

