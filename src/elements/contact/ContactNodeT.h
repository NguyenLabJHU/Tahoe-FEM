/* $Id: ContactNodeT.h,v 1.1 2001-04-16 17:30:50 rjones Exp $ */


#ifndef _CONTACT_NODE_T_H_
#define _CONTACT_NODE_T_H_

/* direct members */
#include "ArrayT.h"
#include "dArrayT.h"
#include "SurfaceT.h"
#include "FaceT.h"

/* forward declarations */
class ofstreamT;

class ContactNodeT 
{
  public:

	/* constructor */
	ContactNodeT(void);

	/* constructor */
	~ContactNodeT(void);

	/* print data */
	void PrintData(ostream& out);

	/* assign opposing point on surface */
	bool AssignOpposing
		(SurfaceT* opposing_surface, FaceT* opposing_face,
		double* xi, double g);
	
  protected:
        /* nodal arrays */
	SurfaceT*  fOpposingSurface ; 
	FaceT*     fOpposingFace ; 
	double     fxi[2] ;
	double     fGap ;

  public:
        /* access functions */ 
        inline const SurfaceT* OpposingSurface(void) {return fOpposingSurface;}
        inline const FaceT* OpposingFace(void) {return fOpposingFace;}
        inline const double* OpposingLocalCoordinates(void) {return fxi;}
        inline const double Gap(void) {return fGap;}


  private:

};

#endif /* _CONTACT_NODE_T_H_ */
