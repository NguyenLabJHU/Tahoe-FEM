/* $Id: ContactSurfaceT.h,v 1.2 2001-04-09 22:28:55 rjones Exp $ */


#ifndef _CONTACT_SURFACE_T_H_
#define _CONTACT_SURFACE_T_H_

/* base class */
#include "SurfaceT.h"

/* direct members */
#include "ArrayT.h"
#include "dArray2DT.h"
#include "dArrayT.h"

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

	/* print data */
	void PrintData(ostream& out);
	
        /* access functions */
        const SurfaceT* OpposingSurface(int node_number) 
		{return fOpposingSurface[node_number];}
        const FaceT* OpposingFace(int node_number) 
		{return fOpposingFace[node_number];}
	double Gap(int node_number)
		{return fGaps[node_number];}


  protected:
        /* nodal arrays */
	ArrayT <SurfaceT*>  fOpposingSurface ; 
	ArrayT <FaceT*>     fOpposingFace ; 
	dArray2DT           fOpposingLocalCoordinates ;
	dArrayT             fGaps ;
  private:

};

#endif /* _CONTACT_SURFACE_T_H_ */
