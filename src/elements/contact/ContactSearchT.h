/* $Id */

#ifndef _CONTACT_SEARCH_T_H_
#define _CONTACT_SEARCH_T_H_

/* base class */
#include "ContactT.h"

/* direct members */
#include "AutoArrayT.h"
#include "nVariArray2DT.h"

/* forward declarations */
class iGridManagerT;

class ContactSearchT: public ContactElementT // ?????
{
public:

	/* constructor */
	ContactSearchT(FEManagerT& fe_manager);

	/* destructor */
	~ContactSearchT(void);

protected:
	/* steps in setting contact configuration */
	void SetInteractions(void); // determine interactions
	void UpdateProjection(void); // updates interactions

private:

	void UpdateKinematicData(void); 

	/* nodes on surface 1 projected onto faces of surface 2*/
	void NodeFaceSearch
		(ContactSurfaceT* surf1, ContactSurfaceT* surf2, parameters); 

	void UpdateProjection
		(ContactSurfaceT* surf1, ContactSurfaceT* surf2, parameters); 
	
protected:
	
	/* search grid */
	iGridManagerT* fGrid;

	MatrixT SearchParameters;

};

#endif /*_CONTACT_SEARCH_T_H_*/
