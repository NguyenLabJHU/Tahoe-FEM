/* $Id: D2MeshFreeSupport2DT.h,v 1.7 2004-07-15 08:30:07 paklein Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D2_MF_SUPPORT_2D_T_H_
#define _D2_MF_SUPPORT_2D_T_H_

/* base class */
#include "D2MeshFreeSupportT.h"

namespace Tahoe {

/** class for support of meshfree field calculations up to second gradients 
 * in two dimensions. See base class documentation for information  about 
 * class initialization. */
 class D2MeshFreeSupport2DT: public D2MeshFreeSupportT
{

public:

	/** constructor.
	 * \param domain used to determine the location of integration points
	 * \param coords array of all particle coordinates 
	 * \param connects integration cell connectivities 
	 * \param nongridnodes index of paricles not included in the connectivities */
	D2MeshFreeSupport2DT(const ParentDomainT* domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes);

	/** set cutting field facets. 
	 * \param facet_coords list of coordinate for each facet: [nfacets] x [num_facet_nodes*nsd] 
	 * \param num_facet_nodes number of nodes defining each facet */
	virtual void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);

private:

	/* process boundaries - nodes marked as "inactive" at the
	 * current x_node by setting dmax = -1.0 */
	virtual void ProcessBoundaries(const dArray2DT& coords,
		const dArrayT& x_node, dArray2DT& nodal_params);

	/* returns 1 if the path x1-x2 is visible */
	virtual int Visible(const double* x1, const double* x2);

	/* returns 1 if the line segment a->b intersects the line
	 * segment p->q */
	int Intersect(const double* a, const double* b, const double* p,
		const double* q) const;
};

} // namespace Tahoe 
#endif /* _D2_MF_SUPPORT_2D_T_H_ */
