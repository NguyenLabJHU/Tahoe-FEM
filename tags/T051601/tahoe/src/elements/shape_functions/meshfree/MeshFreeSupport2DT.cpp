/* $Id: MeshFreeSupport2DT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (09/10/1998)                                          */
/* MLS shape function support for 2D                                      */

#include "MeshFreeSupport2DT.h"

#include <math.h>
#include <string.h>

#include "ExceptionCodes.h"
#include "Constants.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

static    int Max(int a, int b) { return (a > b) ? a : b; };
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
MeshFreeSupport2DT::MeshFreeSupport2DT(const ParentDomainT& domain,
	const dArray2DT& coords, const iArray2DT& connects,
	const iArrayT& nongridnodes, FormulationT code, double dextra,
	int complete, bool store_shape):
	MeshFreeSupportT(domain, coords, connects, nongridnodes, code,
		dextra, complete, store_shape)
{

}

/* cutting facet functions */
void MeshFreeSupport2DT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	/* inherited */
	MeshFreeSupportT::SetCuttingFacets(facet_coords, num_facet_nodes);

	/* checks */
	if (fNumFacetNodes != 2)
	{
		cout << "\n MeshFreeSupport2DT::SetCuttingFacets: 2D cutting facets must\n"
		     <<   "     have 2 nodes: " << fNumFacetNodes << endl;
		throw eSizeMismatch;
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* process boundaries - nodes marked as "inactive" at the
* current x_node by setting dmax = -1.0 */
void MeshFreeSupport2DT::ProcessBoundaries(const dArray2DT& coords,
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
	double eps = 0.01;
	double* pnode = x_node.Pointer();
	for (int j = 0; j < fCutCoords->MajorDim(); j++)
	{
		double* p1 = (*fCutCoords)(j);
		double* p2 = p1 + 2;
		for (int i = 0; i < coords.MajorDim(); i++)
			if (Intersect(pnode, coords(i), p1, p2, eps))
				dmax[i] = -1.0;
	}
}		

/* returns 1 if the path x1-x2 is visible */
int MeshFreeSupport2DT::Visible(const double* x1, const double* x2)
{
	/* quick exit */
	if (!fCutCoords) return 1;

	/* exhaustive search for now */
	double eps = 0.01;
	for (int j = 0; j < fCutCoords->MajorDim(); j++)
	{
		double* p1 = (*fCutCoords)(j);
		double* p2 = p1 + 2;
		if (Intersect(x1, x2, p1, p2, eps)) return 0;
	}
	return 1;
}
	
/* returns 1 if the line segment a->b intersects the line
* segment p->q */
int MeshFreeSupport2DT::Intersect(const double* a, const double* b, const double* p_,
	const double* q_, double eps_pq) const
{
	double xc[2];
	xc[0] = 0.5*(p_[0] + q_[0]);
	xc[1] = 0.5*(p_[1] + q_[1]);
	
	double eps  = eps_pq/2.0;
	double eps1 = 1.0 + eps;
	double p[2];
	p[0] = eps1*p_[0] - eps*xc[0];
	p[1] = eps1*p_[1] - eps*xc[1];
	double q[2];
	q[0] = eps1*q_[0] - eps*xc[0];
	q[1] = eps1*q_[1] - eps*xc[1];

	double v_ab[2] = {b[0] - a[0], b[1] - a[1]};
	double v_pq[2] = {q[0] - p[0], q[1] - p[1]};
	double test;
	
	/* left or right handed intersection */
	test = v_ab[0]*v_pq[1] - v_ab[1]*v_pq[0];
	if (test < 0.0)
	{
		v_ab[0] = -v_ab[0];
		v_ab[1] = -v_ab[1];

		v_pq[0] = -v_pq[0];
		v_pq[1] = -v_pq[1];
	}
	
	double v_ap[2] = {p[0] - a[0], p[1] - a[1]};
	test = v_ap[0]*v_ab[1] - v_ap[1]*v_ab[0];
	if (test < 0.0)
		return 0;
	else
	{
		double v_aq[2] = {q[0] - a[0], q[1] - a[1]}; 	
		test = v_ab[0]*v_aq[1] - v_ab[1]*v_aq[0];
		if (test < 0.0) return 0;
	}
	
	double v_pb[2] = {b[0] - p[0], b[1] - p[1]};
	test = v_pb[0]*v_pq[1] - v_pb[1]*v_pq[0];
	if (test < 0.0)
		return 0;
	else
	{
		double v_pa[2] = {a[0] - p[0], a[1] - p[1]};
		test = v_pq[0]*v_pa[1] - v_pq[1]*v_pa[0];
		if (test < 0.0) return 0;
	}

	return 1;
}

int MeshFreeSupport2DT::Intersect_2(const double* a, const double* b, const double* p_,
	const double* q_, double eps_pq) const
{
	/* extend segment pq by 2*eps */
	double xc[2];
	xc[0] = 0.5*(p_[0] + q_[0]);
	xc[1] = 0.5*(p_[1] + q_[1]);
	
	double eps  = eps_pq/2.0;
	double eps1 = 1.0 + eps;
	double p[2];
	p[0] = eps1*p_[0] - eps*xc[0];
	p[1] = eps1*p_[1] - eps*xc[1];
	double q[2];
	q[0] = eps1*q_[0] - eps*xc[0];
	q[1] = eps1*q_[1] - eps*xc[1];

	/* load */
	double a1 = a[0];
	double a2 = a[1];
	double b1 = b[0];
	double b2 = b[1];

	double p1 = p[0];
	double p2 = p[1];
	double q1 = q[0];
	double q2 = q[1];

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20;

	z1 = -a2*b1;
	z2 = a1*b2;
	z3 = -a2*p1;
	z4 = a2*p1;
	z5 = -b2*p1;
	z6 = b2*p1;
	z7 = -a1*p2;
	z8 = a1*p2;
	z9 = -b1*p2;
	z10 = b1*p2;
	z11 = -a2*q1;
	z12 = a2*q1;
	z13 = -b2*q1;
	z14 = b2*q1;
	z15 = -p2*q1;
	z16 = -a1*q2;
	z17 = a1*q2;
	z18 = -b1*q2;
	z19 = b1*q2;
	z20 = p1*q2;
	z4 = z10 + z4 + z5 + z7;
	z5 = z11 + z14 + z17 + z18 + z4;
	z1 = z1 + z2 + z4;
	z2 = z12 + z15 + z16 + z20 + z3 + z8;
	z3 = z12 + z13 + z16 + z19 + z3 + z6 + z8 + z9;
	
	/* return - {s_ab, s_pq}
	 *
	 *   s_ab = z2/z3
	 *   s_pq = z1/z5
	 *
	 * intersect if both are [0 1] */

	//z4 = 1./z5;
	//z3 = 1./z3;
	//z1 = z1*z4;
	//z2 = z2*z3;
	//List(z2,z1);
	
	if (fabs(z3) < kSmall || fabs(z4) < kSmall)
		return 0;
	else
	{
		double s_ab = z2/z3;
		double s_pq = z1/z5;
		if (s_ab < 0 || s_ab > 1 ||
		    s_pq < 0 || s_pq > 1)
		    return 0;
		else
			return 1;
	}
}
