/* $Id: Contact3DT.h,v 1.5 2003-03-02 18:57:31 paklein Exp $ */
/* created: paklein (07/17/1999) */
#ifndef _CONTACT3D_T_H_
#define _CONTACT3D_T_H_

/* base class */
#include "ContactT.h"

/* direct members */
#include "AutoArrayT.h"
#include "nVariArray2DT.h"

namespace Tahoe {

/* forward declarations */
class iGridManager3DT;

/** base class for 3D contact elements */
class Contact3DT: public ContactT
{
public:

	/* constructor */
	Contact3DT(const ElementSupportT& support, const FieldT& field);

	/* destructor */
	virtual ~Contact3DT(void);

protected:

	/* element data */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/** \name steps in setting contact configuration */
	/*@{*/
	/** set "internal" data */
	virtual bool SetActiveInteractions(void);

	/** set "external" data - interface to FEManager */
	virtual void SetConnectivities(void);
	/*@}*/

	/* convert quad facets to tri's */
	void ConvertQuadToTri(iArray2DT& surface) const;

	/** set surface normal derivative matrix */
	void Set_dn_du(const dArray2DT& curr_coords, dMatrixT& dn_du) const;

	/** second variation of gap vector for 3-noded facet */
	void DDg_tri_facet(
		double* X1, double* X2, double* X3, double* pX,
		double* d1, double* d2, double* d3, double* pd,
		dMatrixT& K) const;

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

} // namespace Tahoe 
#endif /* _CONTACT3D_T_H_ */
