/* $Id: Contact3DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (07/17/1999)                                          */

#ifndef _CONTACT3D_T_H_
#define _CONTACT3D_T_H_

/* base class */
#include "ContactT.h"

/* direct members */
#include "AutoArrayT.h"
#include "nVariArray2DT.h"

/* forward declarations */
class iGridManager3DT;

class Contact3DT: public ContactT
{
public:

	/* constructor */
	Contact3DT(FEManagerT& fe_manager);

	/* destructor */
	virtual ~Contact3DT(void);

protected:

	/* element data */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/* steps in setting contact configuration */
	virtual bool SetActiveInteractions(void); // "internal" data
	virtual void SetConnectivities(void); // "external" data - interface to FEManager

	/* convert quad facets to tri's */
	void ConvertQuadToTri(iArray2DT& surface) const;

private:

	/* sets active striker data (based on current bodies data) */
	void SetActiveStrikers(void); // one contact per striker

	/* return true if vector from A to B intersects the facet */
	bool Intersect(const dArrayT& x1, const dArrayT& x2, const dArrayT& x3,
		const dArrayT& xs, double& h) const;
	
protected:
	
	/* search grid */
	iGridManager3DT* fGrid3D;

	/* work space */
	dArrayT	fx1, fx2, fx3; // facet node coords (shallow)
	dArrayT fStriker;      // striker node coords (shallow)
};

#endif /* _CONTACT3D_T_H_ */
