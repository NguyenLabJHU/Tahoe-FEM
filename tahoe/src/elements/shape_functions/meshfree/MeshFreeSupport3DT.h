/* $Id: MeshFreeSupport3DT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (09/13/1998)                                          */
/* meshfree shape function support for 3D                                 */

#ifndef _MF_SUPPORT_3D_T_H_
#define _MF_SUPPORT_3D_T_H_

/* base class */
#include "MeshFreeSupportT.h"

class MeshFreeSupport3DT: public MeshFreeSupportT
{
public:

	/* constructor */
	MeshFreeSupport3DT(const ParentDomainT& domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes, FormulationT code,
		double dextra, int complete, bool store_shape);

	/* cutting facet functions */
	virtual void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);

private:

	/* process boundaries - nodes marked as "inactive" at the
	 * current x_node by setting dmax = -1.0 */
	virtual void ProcessBoundaries(const dArray2DT& coords,
		const dArrayT& x_node, dArrayT& dmax);

	/* returns 1 if the path x1-x2 is visible */
	virtual int Visible(const double* x1, const double* x2);
	
	/* returns 1 if the line segment a->b intersects the triangular facet
	 * defined by {x0, x1, x2} */
	int Intersect(const double* x0, const double* x1, const double* x2,
		const double *vA, const double* vB) const;

	/* TEMP: fixed cutting surface shapes */
	void CutCircle(const dArray2DT& coords, const dArrayT& x_node,
		dArrayT& dmax);

	void CutEllipse(const dArray2DT& coords, const dArrayT& x_node,
		dArrayT& dmax);
};

#endif /* _MF_SUPPORT_3D_T_H_ */
