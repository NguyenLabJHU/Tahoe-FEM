/* $Id: ContactSearchT.h,v 1.11 2003-11-21 22:54:34 paklein Exp $ */
#ifndef _CONTACT_SEARCH_T_H_
#define _CONTACT_SEARCH_T_H_

/* direct members */
#include "nMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iGridManagerT.h"

namespace Tahoe {

/* forward declarations */
class FEManagerT;
class ContactSurfaceT;

class ContactSearchT
{
public:

	/* constructor */
	ContactSearchT(ArrayT<ContactSurfaceT>& surfaces,
		nMatrixT<dArrayT>& search_parameters);

	/* destructor */
	~ContactSearchT(void);

	/* determines contact configuration */
	bool SetInteractions(void);

	/* updates contact configuration */
	bool UpdateInteractions(void); 
protected:

private:
	/* does initization and tracking */
	void Initialize(void); 

	/* nodes on surface 1 projected onto faces of surface 2*/
	void NodeFaceSearch
		(ContactSurfaceT& surface1, ContactSurfaceT& surface2,
		const dArrayT& parameters); 

	/* update gaps and local coordinates of projection */
	bool UpdateProjection(void);
	
	/* search grid */
	iGridManagerT* fGrid;

	/* surface (data) */
	ArrayT<ContactSurfaceT>& fSurfaces;

	/* search parameters from contact element */
	const nMatrixT<dArrayT>& fSearchParameters;

	/*workspace*/
	double centroid[3];
	double radius;
	iArrayT grid_nodes;

};

} // namespace Tahoe 
#endif /*_CONTACT_SEARCH_T_H_*/
