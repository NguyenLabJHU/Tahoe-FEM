/* $Id: ContactSearchT.h,v 1.2 2001-04-09 22:28:55 rjones Exp $ */

#ifndef _CONTACT_SEARCH_T_H_
#define _CONTACT_SEARCH_T_H_


/* direct members */
#include "nMatrixT.h"
#include "dArrayT.h"
#include "iGridManagerT.h"

/* forward declarations */
class FEManagerT;
class ContactSurfaceT;

class ContactSearchT
{
public:

	/* constructor */
	ContactSearchT(FEManagerT& fe_manager, 
		ArrayT<ContactSurfaceT>& surfaces,
		nMatrixT<dArrayT>& search_parameters);

	/* destructor */
	~ContactSearchT(void);

	/* determines contact configuration */
	bool SetInteractions(void);

	/* updates contact configuration */
	bool UpdateInteractions(void); 
protected:

private:

	void Initialize(void); 

	/* nodes on surface 1 projected onto faces of surface 2*/
	void NodeFaceSearch
		(ContactSurfaceT& surface1, ContactSurfaceT& surface2); 

	/* update gaps and local coordinates of projection */
	void UpdateProjection
		(ContactSurfaceT& surface1, ContactSurfaceT& surface2); 
	
	/* search grid */
	iGridManagerT* fGrid;

	/* surface (data) */
	const ArrayT<ContactSurfaceT>& fSurfaces;

	/* search parameters from contact element */
	const nMatrixT<dArrayT>& fSearchParameters;

};

#endif /*_CONTACT_SEARCH_T_H_*/
