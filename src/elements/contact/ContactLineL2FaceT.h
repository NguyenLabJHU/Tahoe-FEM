/* $Id: ContactLineL2FaceT.h,v 1.2 2002-07-02 19:55:19 cjkimme Exp $ */

#ifndef _CONTACT_LINEL2_FACE_T_H_
#define _CONTACT_LINEL2_FACE_T_H_

/* base class */
#include "ContactFaceT.h"


namespace Tahoe {

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

} // namespace Tahoe 
#endif /* _CONTACT_LINEL2_FACE_T_H_ */

