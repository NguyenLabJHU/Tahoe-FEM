/* $Id: LatLongPtsT.h,v 1.3.56.1 2004-06-19 23:28:07 paklein Exp $ */
/* created: paklein (10/31/1997) */
#ifndef _LATLONG_PTS_T_H_
#define _LATLONG_PTS_T_H_

/* base class */
#include "SpherePointsT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

class LatLongPtsT: public SpherePointsT
{
public:

	/** constructor */
	LatLongPtsT(int n_phi, int n_theta);

	/** generate sphere points:
	 *
	 *   phi   = angle about z from x
	 *   theta = angle about x from z
	 *
	 * The final orientation is generated by applied the
	 * phi and theta rotations in succession about the local
	 * axes.
	 */
	virtual const dArray2DT& SpherePoints(double phi, double theta);

private:

	/* parameters */
	int	fNphi;
	int fNtheta;
};

} // namespace Tahoe 
#endif /* _LATLONG_PTS_T_H_ */
