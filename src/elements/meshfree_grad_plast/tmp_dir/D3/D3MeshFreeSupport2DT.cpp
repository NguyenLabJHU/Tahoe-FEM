/* $Id: D3MeshFreeSupport2DT.cpp,v 1.1 2004-09-28 23:13:59 raregue Exp $ */
/* created: paklein (10/23/1999) */
#include "D3MeshFreeSupport2DT.h"

#include <math.h>
#include <string.h>

#include "ExceptionT.h"
#include "toolboxConstants.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

using namespace Tahoe;

static    int Max(int a, int b) { return (a > b) ? a : b; };
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
D3MeshFreeSupport2DT::D3MeshFreeSupport2DT(const ParentDomainT* domain,
	const dArray2DT& coords, const iArray2DT& connects, const iArrayT& nongridnodes):
	D3MeshFreeSupportT(domain, coords, connects, nongridnodes)
{

}

/* cutting facet functions */
void D3MeshFreeSupport2DT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	/* inherited */
	D3MeshFreeSupportT::SetCuttingFacets(facet_coords, num_facet_nodes);

	/* checks */
	if (fNumFacetNodes != 2)
		ExceptionT::SizeMismatch("D3MeshFreeSupport2DT::SetCuttingFacets", "2D cutting facets must have 2 nodes: %d",
			fNumFacetNodes);
}

/*************************************************************************
* Private
*************************************************************************/

/* process boundaries - nodes marked as "inactive" at the
* current x_node by setting dmax = -1.0 */
void D3MeshFreeSupport2DT::ProcessBoundaries(const dArray2DT& coords,
	const dArrayT& x_node, dArray2DT& nodal_params)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (coords.MajorDim() != nodal_params.MajorDim()) throw ExceptionT::kSizeMismatch;
	if (coords.MinorDim() != x_node.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* quick exit */
	if (!fCutCoords) return;
	
	/* exhaustive search for now */
	const double* pnode = x_node.Pointer();
	for (int j = 0; j < fCutCoords->MajorDim(); j++)
	{
		const double* p1 = (*fCutCoords)(j);
		const double* p2 = p1 + 2;
		for (int i = 0; i < coords.MajorDim(); i++)
			if (Intersect(p1, p2, pnode, coords(i))) 
				nodal_params.SetRow(i, -1.0);
	}
}		

/* returns 1 if the path x1-x2 is visible */
int D3MeshFreeSupport2DT::Visible(const double* x1, const double* x2)
{
	/* quick exit */
	if (!fCutCoords) return 1;

	/* exhaustive search for now */
	for (int j = 0; j < fCutCoords->MajorDim(); j++)
	{
		const double* p1 = (*fCutCoords)(j);
		const double* p2 = p1 + 2;
		if (Intersect(x1, x2, p1, p2)) return 0;
	}
	return 1;
}
	
/* returns 1 if the line segment a->b intersects the line
* segment p->q */
int D3MeshFreeSupport2DT::Intersect(const double* a, const double* b, const double* p,
	const double* q) const
{
	double v_ab[2] = {b[0] - a[0], b[1] - a[1]};
	double v_pq[2] = {q[0] - p[0], q[1] - p[1]};
	
	/* left or right handed intersection */
	if (v_ab[0]*v_pq[1] - v_ab[1]*v_pq[0] < kSmall)
	{
		v_ab[0] = -v_ab[0];
		v_ab[1] = -v_ab[1];

		v_pq[0] = -v_pq[0];
		v_pq[1] = -v_pq[1];
	}
	
	double v_ap[2] = {p[0] - a[0], p[1] - a[1]};
	if (v_ap[0]*v_ab[1] - v_ap[1]*v_ab[0] < kSmall)
		return 0;
	else
	{
		double v_aq[2] = {q[0] - a[0], q[1] - a[1]}; 	
		if (v_ab[0]*v_aq[1] - v_ab[1]*v_aq[0] < kSmall) return 0;
	}
	
	double v_pb[2] = {b[0] - p[0], b[1] - p[1]};
	if (v_pb[0]*v_pq[1] - v_pb[1]*v_pq[0] < kSmall)
		return 0;
	else
	{
		double v_pa[2] = {a[0] - p[0], a[1] - p[1]};
		if ( v_pq[0]*v_pa[1] - v_pq[1]*v_pa[0] < kSmall ) return 0;
	}

	return 1;
}
