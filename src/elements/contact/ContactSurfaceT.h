/*
 * File: ContactSurfaceT.h
 */


#ifndef _CONTACT_SURFACE_T_H_
#define _CONTACT_SURFACE_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class FEManagerT;

/* direct members */
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

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

	/* print data */
	void PrintData(ostream& out);

  protected:
        /* nodal arrays */
	ArrayT <SurfaceT*>  fOpposingSurface ; 
	ArrayT <FaceT*>     fOpposingFace ; 
	dArray2DT           fOpposingLocalCoordinates ;
	dArrayT             fGaps ;
  private:

};

#endif /* _CONTACT_SURFACE_T_H_ */
