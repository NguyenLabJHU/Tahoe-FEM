/* $Id: MeshFreeSupport2DT.h,v 1.8 2004-07-15 08:29:59 paklein Exp $ */
/* created: paklein (09/10/1998) */

#ifndef _MF_SUPPORT_2D_T_H_
#define _MF_SUPPORT_2D_T_H_

/* base class */
#include "MeshFreeSupportT.h"

namespace Tahoe {

/** Class for support of meshfree methods in two dimensions. See
 * documentation from base class for information about initialization. */
class MeshFreeSupport2DT: public MeshFreeSupportT
{
public:

	/** constructor.
	 * \param domain used to determine the location of integration points
	 * \param coords array of all particle coordinates 
	 * \param connects integration cell connectivities 
	 * \param nongridnodes index of paricles not included in the connectivities
	 * \param in input stream for class and window function parameters */
	MeshFreeSupport2DT(const ParentDomainT* domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes);

	/** construct object sufficient for calling methods inherited from ParameterInterfaceT
	 * to collect the class parameters, but not for doing any meshfree calculations */
	MeshFreeSupport2DT(void);

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

	/* returns 1 if the line segment a->b intersects the line
	 * segment p->q, extending pq by eps to create overlap */
	int Intersect(const double* a, const double* b, const double* p,
		const double* q, double eps_pq) const;

	/* different approach - extending pq by eps to create overlap*/
	int Intersect_2(const double* a, const double* b, const double* p,
		const double* q, double eps_pq) const;
};

} // namespace Tahoe 
#endif /* _MF_SUPPORT_2D_T_H_ */
