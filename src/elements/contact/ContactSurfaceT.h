/* $Id: ContactSurfaceT.h,v 1.8 2001-06-27 18:16:21 rjones Exp $ */


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

	/* potential connectivities based on growing/sliding contact */
	void SetPotentialConnectivity(void);

	/* access functions */
	inline ArrayT<ContactNodeT*>& ContactNodes(void) 
		{return fContactNodes;}
	inline RaggedArray2DT<int>& Connectivities(void)
		{return fConnectivities;}
	inline RaggedArray2DT<int>& EqNums(void)
		{return fEqNums;}  // this can NOT be const
#if 0
	inline ArrayT<ContactNodeT*>& PreviousContactNodes(void) 
		{return fPreviousContactNodes;}
#endif
#if 0
	void PrintContactArea(ofstream& out) const;
#endif
	void PrintContactArea(ostream& out) const;
	void PrintGap(ostream& out) const;
	void PrintGap(ofstream& out) const;
	void PrintNormals(ofstream& out) const;

  protected:
        /* nodal arrays */
	ArrayT <ContactNodeT*>  fContactNodes ; 

	/* potential connectivities for the time step */
	RaggedArray2DT<int> fConnectivities;

	/* space for associated equation numbers */
	RaggedArray2DT<int> fEqNums;
#if 0
	/* for frictional slip */
	ArrayT <ContactNodeT*>  fPreviousContactNodes;
#endif

};

#endif /* _CONTACT_SURFACE_T_H_ */
