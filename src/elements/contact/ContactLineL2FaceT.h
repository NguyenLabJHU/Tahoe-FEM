/* $Id: ContactLineL2FaceT.h,v 1.1 2001-09-24 20:43:24 rjones Exp $ */

#ifndef _CONTACT_LINEL2_FACE_T_H_
#define _CONTACT_LINEL2_FACE_T_H_

/* base class */
#include "ContactFaceT.h"

class ContactLineL2FaceT : public ContactFaceT
{
public:

	/* constructor */
	ContactLineL2FaceT (FaceT* face); 	

	/* destructor */
  	~ContactLineL2FaceT(void);

	void ComputePressureFunctions
		(const double* local_coordinates, dMatrixT& shape_functions)
        const;

};

#endif /* _CONTACT_LINEL2_FACE_T_H_ */

