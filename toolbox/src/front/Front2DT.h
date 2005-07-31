/* $Id: Front2DT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/18/1999)                                          */

#ifndef _FRONT_2D_T_H_
#define _FRONT_2D_T_H_

/* base class */
#include "FrontT.h"

/* forward declarations */
class FrontSegmentT;

class Front2DT: public FrontT
{
public:

	/* constructor */
	Front2DT(double cone, double da, double da_s, int num_pts);

	/* construct initial front */
	virtual void Initialize(const dArray2DT& facet_coords, const iArrayT& fr_facets,
		const iArrayT& fr_edges);
	
	/* extend front at the specified points along the given direction - returns new
	 * cutting facets */
	virtual const dArray2DT& NewFacets(const ArrayT<int>& extend_pts,
		const ArrayT<int>& extend_dir);
};

#endif /* _FRONT_2D_T_H_ */
