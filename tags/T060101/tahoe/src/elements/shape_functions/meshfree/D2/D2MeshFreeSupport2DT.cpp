/* $Id: D2MeshFreeSupport2DT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/23/1999)                                          */

#include "D2MeshFreeSupport2DT.h"

#include <math.h>
#include <string.h>

#include "ExceptionCodes.h"
#include "Constants.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

static    int Max(int a, int b) { return (a > b) ? a : b; };
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
D2MeshFreeSupport2DT::D2MeshFreeSupport2DT(const ParentDomainT& domain,
	const dArray2DT& coords, const iArray2DT& connects, const iArrayT& nongridnodes,
	FormulationT code, double dextra, int complete, bool store_shape):
	D2MeshFreeSupportT(domain, coords, connects, nongridnodes, code, dextra,
		complete, store_shape)
{

}

/* cutting facet functions */
void D2MeshFreeSupport2DT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	/* inherited */
	D2MeshFreeSupportT::SetCuttingFacets(facet_coords, num_facet_nodes);

	/* checks */
	if (fNumFacetNodes != 2)
	{
		cout << "\n D2MeshFreeSupport2DT::SetCuttingFacets: 2D cutting facets must\n"
		     <<   "     have 2 nodes: " << fNumFacetNodes << endl;
		throw eSizeMismatch;
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* process boundaries - nodes marked as "inactive" at the
* current x_node by setting dmax = -1.0 */
void D2MeshFreeSupport2DT::ProcessBoundaries(const dArray2DT& coords,
	const dArrayT& x_node, dArrayT& dmax)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (coords.MajorDim() != dmax.Length()) throw eSizeMismatch;
	if (coords.MinorDim() != x_node.Length()) throw eSizeMismatch;
#endif

	/* quick exit */
	if (!fCutCoords) return;
	
	/* exhaustive search for now */
	double* pnode = x_node.Pointer();
	for (int j = 0; j < fCutCoords->MajorDim(); j++)
	{
		double* p1 = (*fCutCoords)(j);
		double* p2 = p1 + 2;
	
		for (int i = 0; i < coords.MajorDim(); i++)
			if (Intersect(p1, p2, pnode, coords(i))) dmax[i] = -1.0;
	}
}		

/* returns 1 if the path x1-x2 is visible */
int D2MeshFreeSupport2DT::Visible(const double* x1, const double* x2)
{
	/* quick exit */
	if (!fCutCoords) return 1;

	/* exhaustive search for now */
	for (int j = 0; j < fCutCoords->MajorDim(); j++)
	{
		double* p1 = (*fCutCoords)(j);
		double* p2 = p1 + 2;
	
		if (Intersect(x1, x2, p1, p2)) return 0;
	}
	return 1;
}
	
/* returns 1 if the line segment a->b intersects the line
* segment p->q */
int D2MeshFreeSupport2DT::Intersect(const double* a, const double* b, const double* p,
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
