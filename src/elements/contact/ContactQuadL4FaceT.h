/* $Id: ContactQuadL4FaceT.h,v 1.1 2001-09-24 20:43:25 rjones Exp $ */

#ifndef _CONTACT_QUADL4_FACE_T_H_
#define _CONTACT_QUADL4_FACE_T_H_

/* base class */
#include "ContactFaceT.h"

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

#endif /* _CONTACT_QUADL4_FACE_T_H_ */

