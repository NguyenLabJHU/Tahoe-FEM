/* $Id: ContactLineQ3FaceT.h,v 1.1 2001-09-24 20:43:25 rjones Exp $ */

#ifndef _CONTACT_LINEQ3_FACE_T_H_
#define _CONTACT_LINEQ3_FACE_T_H_

/* base class */
#include "ContactFaceT.h"

class ContactLineQ3FaceT : public ContactFaceT
{
public:

	/* constructor */
	ContactLineQ3FaceT (FaceT* face); 	

	/* destructor */
  	~ContactLineQ3FaceT(void);

	void ComputePressureFunctions
		(const double* local_coordinates, dMatrixT& shape_functions)
        const;

};

#endif /* _CONTACT_LINEQ3_FACE_T_H_ */

