/* $Id: ContactSurfaceT.h,v 1.4 2001-04-16 17:30:51 rjones Exp $ */


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

	/* constructor */
	~ContactSurfaceT(void);

	/* allocate contact node array */
	void AllocateContactNodes(void);

  protected:
        /* nodal arrays */
	ArrayT <ContactNodeT>  fContactNodes ; 

};

#endif /* _CONTACT_SURFACE_T_H_ */
