/* $Id: ContactSurfaceT.h,v 1.5 2001-04-19 23:47:01 rjones Exp $ */


#ifndef _CONTACT_SURFACE_T_H_
#define _CONTACT_SURFACE_T_H_

/* base class */
#include "SurfaceT.h"

/* direct members */
#include "ArrayT.h"
#include "SurfaceT.h"
#include "ContactNodeT.h"

/* forward declarations */
class ofstreamT;
class FEManagerT;

/* 
a ContactSurface will only have one opposing face per
node and be considered "smooth" i.e. the full boundary 
surface of a cube will be made up of 6 surfaces
*/

class ContactSurfaceT : public SurfaceT
{
  public:
	/* constructor */
	ContactSurfaceT(void);

	/* destructor */
	~ContactSurfaceT(void);

	/* allocate contact node array */
	void AllocateContactNodes(void);

	/* move current to previous */
	void CopyCurrentToPrevious(void);

	/* access functions */
	inline ArrayT<ContactNodeT*>& ContactNodes(void) 
		{return fContactNodes;}
#if 0
	inline ArrayT<ContactNodeT*>& PreviousContactNodes(void) 
		{return fPreviousContactNodes;}
#endif


  protected:
        /* nodal arrays */
	ArrayT <ContactNodeT*>  fContactNodes ; 

#if 0
	/* for frictional slip */
	ArrayT <ContactNodeT*>  fPreviousContactNodes;
#endif

};

#endif /* _CONTACT_SURFACE_T_H_ */
