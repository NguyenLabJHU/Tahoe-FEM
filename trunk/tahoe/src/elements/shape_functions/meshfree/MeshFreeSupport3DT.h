/* $Id: MeshFreeSupport3DT.h,v 1.7 2004-01-27 01:21:11 cjkimme Exp $ */
/* created: paklein (09/13/1998) */

#ifndef _MF_SUPPORT_3D_T_H_
#define _MF_SUPPORT_3D_T_H_

/* base class */
#include "MeshFreeSupportT.h"

namespace Tahoe {

/** Class for support of meshfree methods in three dimensions. See
 * documentation from base class for information about initialization. */
class MeshFreeSupport3DT: public MeshFreeSupportT
{
public:

	/** constructor.
	 * \param domain used to determine the location of integration points
	 * \param coords array of all particle coordinates 
	 * \param connects integration cell connectivities 
	 * \param nongridnodes index of paricles not included in the connectivities
	 * \param in input stream for class and window function parameters */
	MeshFreeSupport3DT(const ParentDomainT* domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes, ifstreamT& in);

	/** set field cutting facets. 
	 * \param facet_coords list of coordinate for each facet: [nfacets] x [num_facet_nodes*nsd] 
	 * \param num_facet_nodes number of nodes defining each facet */
	virtual void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);

private:

	/* process boundaries - nodes marked as "inactive" at the
	 * current x_node by setting nodal_params = -1.0 */
	virtual void ProcessBoundaries(const dArray2DT& coords,
		const dArrayT& x_node, dArray2DT& nodal_params);

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

} // namespace Tahoe 
#endif /* _MF_SUPPORT_3D_T_H_ */
