/* $Id: ContactLineQ3FaceT.h,v 1.2 2002-07-02 19:55:19 cjkimme Exp $ */

#ifndef _CONTACT_LINEQ3_FACE_T_H_
#define _CONTACT_LINEQ3_FACE_T_H_

/* base class */
#include "ContactFaceT.h"


namespace Tahoe {

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

} // namespace Tahoe 
#endif /* _CONTACT_LINEQ3_FACE_T_H_ */

