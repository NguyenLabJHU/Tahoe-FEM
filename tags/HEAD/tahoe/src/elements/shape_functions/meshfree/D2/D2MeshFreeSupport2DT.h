/* $Id: D2MeshFreeSupport2DT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/23/1999)                                          */

#ifndef _D2_MF_SUPPORT_2D_T_H_
#define _D2_MF_SUPPORT_2D_T_H_

/* base class */
#include "D2MeshFreeSupportT.h"

class D2MeshFreeSupport2DT: public D2MeshFreeSupportT
{
public:

	/* constructor */
	D2MeshFreeSupport2DT(const ParentDomainT& domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes, FormulationT code,
		double dextra, int complete, bool storeshape);

	/* cutting facet functions */
	virtual void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);

private:

	/* process boundaries - nodes marked as "inactive" at the
	 * current x_node by setting dmax = -1.0 */
	virtual void ProcessBoundaries(const dArray2DT& coords,
		const dArrayT& x_node, dArrayT& dmax);

	/* returns 1 if the path x1-x2 is visible */
	virtual int Visible(const double* x1, const double* x2);

	/* returns 1 if the line segment a->b intersects the line
	 * segment p->q */
	int Intersect(const double* a, const double* b, const double* p,
		const double* q) const;
};

#endif /* _D2_MF_SUPPORT_2D_T_H_ */
