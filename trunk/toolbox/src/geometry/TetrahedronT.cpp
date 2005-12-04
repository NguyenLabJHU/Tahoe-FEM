/* $Id: TetrahedronT.cpp,v 1.10 2005-12-04 16:56:29 paklein Exp $ */
/* created: paklein (10/22/1996) */
#include "TetrahedronT.h"
#include "QuadT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* parameters */
const int kTetnsd           = 3;
const int kNumVertexNodes	= 4;
const int kNumFacets        = 4;

/* constructor */
TetrahedronT::TetrahedronT(int numnodes): GeometryBaseT(numnodes, kNumFacets) {}

/* evaluate the shape functions and gradients. */
void TetrahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
	const char caller[] = "TetrahedronT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
	if (coords.Length() != 3 ||
	        Na.Length() != fNumNodes) ExceptionT::SizeMismatch(caller);
	if (fNumNodes != 4 && fNumNodes != 10) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
#if 0
	if (r < 0.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < 0.0 || s > 1.0) ExceptionT::OutOfRange(caller);
	if (t < 0.0 || t > 1.0) ExceptionT::OutOfRange(caller);
#endif

	if (fNumNodes == 4)
	{
		/* shape functions */
		Na[0] = r;
		Na[1] = s;
		Na[3] = t;
		Na[2] = 1.0 - r - s - t;
	}
	else if (fNumNodes == 10)
	{
		double u = 1.0 - r - s - t;
	
		/* shape functions */
		Na[0] = r*(2.0*r - 1.0);
		Na[1] = s*(2.0*s - 1.0);
		Na[2] = u*(2.0*u - 1.0);
		Na[3] = t*(2.0*t - 1.0);

		Na[4] = 4.0*r*s;
		Na[5] = 4.0*s*u;
		Na[6] = 4.0*r*u;
		Na[7] = 4.0*r*t;
		Na[8] = 4.0*s*t;
		Na[9] = 4.0*t*u;
	}
	else
		ExceptionT::GeneralFail(caller);
}

/* evaluate the shape functions and gradients. */
void TetrahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
{
	const char caller[] = "TetrahedronT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
	if (coords.Length() != 3 ||
	        Na.Length() != fNumNodes ||
	     DNa.MajorDim() != 3 ||
	     DNa.MinorDim() != fNumNodes) ExceptionT::SizeMismatch(caller);
	if (fNumNodes != 4 && fNumNodes != 10) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
#if 0
	if (r < 0.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < 0.0 || s > 1.0) ExceptionT::OutOfRange(caller);
	if (t < 0.0 || t > 1.0) ExceptionT::OutOfRange(caller);
#endif

	if (fNumNodes == 4)
	{
		/* shape functions */
		Na[0] = r;
		Na[1] = s;
		Na[3] = t;
		Na[2] = 1.0 - r - s - t;

		/* derivatives */
		double* nax = DNa(0);
		double* nay = DNa(1);
		double* naz = DNa(2);

		/* Na,r */
		nax[0] = 1.0;
		nax[1] = 0.0;
		nax[3] = 0.0;
		nax[2] =-1.0;
	
		/* Na,s */
		nay[0] = 0.0;
		nay[1] = 1.0;
		nay[3] = 0.0;
		nay[2] =-1.0;

		/* Na,t */
		naz[0] = 0.0;
		naz[1] = 0.0;
		naz[3] = 1.0;
		naz[2] =-1.0;
	}
	else if (fNumNodes == 10)
	{
		double u = 1.0 - r - s - t;
	
		/* shape functions */
		Na[0] = r*(2.0*r - 1.0);
		Na[1] = s*(2.0*s - 1.0);
		Na[2] = u*(2.0*u - 1.0);
		Na[3] = t*(2.0*t - 1.0);

		Na[4] = 4.0*r*s;
		Na[5] = 4.0*s*u;
		Na[6] = 4.0*r*u;
		Na[7] = 4.0*r*t;
		Na[8] = 4.0*s*t;
		Na[9] = 4.0*t*u;

		/* derivatives */
		double* nax = DNa(0);
		double* nay = DNa(1);
		double* naz = DNa(2);

		/* Na,r */
		nax[0] = 4.0*r - 1.0;
		nax[1] = 0.0;
		nax[2] = 4.0*(r + s + t) - 3.0;
		nax[3] = 0.0;

		nax[4] = 4.0*s;
		nax[5] =-4.0*s;
		nax[6] =-4.0*(2.0*r + s + t - 1.0);
		nax[7] = 4.0*t;
		nax[8] = 0.0;
		nax[9] =-4.0*t;
	
		/* Na,s */
		nay[0] = 0.0;
		nay[1] = 4.0*s - 1.0;
		nay[2] = 4.0*(r + s + t) - 3.0;
		nay[3] = 0.0;

		nay[4] = 4.0*r;
		nay[5] =-4.0*(r + 2.0*s + t - 1.0);
		nay[6] =-4.0*r;
		nay[7] = 0.0;
		nay[8] = 4.0*t;
		nay[9] =-4.0*t;

		/* Na,t */
		naz[0] = 0.0;
		naz[1] = 0.0;
		naz[2] = 4.0*(r + s + t) - 3.0;
		naz[3] = 4.0*t - 1.0;

		naz[4] = 0.0;
		naz[5] =-4.0*s;
		naz[6] =-4.0*r;
		naz[7] = 4.0*r;
		naz[8] = 4.0*s;
		naz[9] =-4.0*(r + s + 2.0*t - 1.0);
	}
	else
		ExceptionT::GeneralFail(caller);
}	

/* compute local shape functions and derivatives */
void TetrahedronT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
	const char caller[] = "TetrahedronT::SetLocalShape";

	/* dimensions */
	int numnodes = Na.MinorDim();
	int numint   = weights.Length();
	int nsd      = Na_x[0].MajorDim();

	/* dimension checks */
	if (numnodes != 4 && numnodes != 10)
		ExceptionT::GeneralFail(caller, "unsupported number of element nodes: %d", numnodes);
	
	if (numint != 1 && numint != 4)
		ExceptionT::GeneralFail(caller, "unsupported number of integration points: %d", numint);

	if (nsd != kTetnsd) ExceptionT::GeneralFail(caller);

	/* initialize */
	Na = 0.0;
	for (int j = 0; j < Na_x.Length(); j++)
		Na_x[j] = 0.0;

	/* integration point coordinates */
	
	/* 1 point */
	double r1[1] = {0.25};
	double s1[1] = {0.25};
	double t1[1] = {0.25};
	
	/* 4 point */
	double r4[4] = {0.58541020, 0.13819660, 0.13819660, 0.13819660};
	double s4[4] = {0.13819660, 0.58541020, 0.13819660, 0.13819660};
	double t4[4] = {0.13819660, 0.13819660, 0.13819660, 0.58541020};
		
	double* r;
	double* s;
	double* t;

	/* set weights and r,s,t */
	switch (numint)
	{
		case 1:	
		{	
			weights[0] = 1.0/6.0;

			/* set coordinates */
			r = r1;
			s = s1;
			t = t1;
		
			break;
		}
		case 4:
		{
			weights = 1.0/24.0;
		
			/* set coordinates */
			r = r4;
			s = s4;
			t = t4;
			
			break;
		}
		default:
			ExceptionT::GeneralFail(caller);
	}	

	/* shape functions and derivatives */
	dArrayT Na_i;
	dArrayT ip_coords(3);
	for (int i = 0; i < numint; i++)
	{
		/* integration point coordinates */
		ip_coords[0] = r[i];
		ip_coords[1] = s[i];
		ip_coords[2] = t[i];
	
		/* compute shape functions and derivatives */
		Na.RowAlias(i, Na_i);
		TetrahedronT::EvaluateShapeFunctions(ip_coords, Na_i, Na_x[i]);
	}
}

/* set the values of the nodal extrapolation matrix */
void TetrahedronT::SetExtrapolation(dMatrixT& extrap) const
{
	const char caller[] = "TetrahedronT::SetExtrapolation";

	/* dimensions */
	int numnodes = extrap.Rows();
	int numint   = extrap.Cols();

	/* dimension checks */
	if (numnodes != 4 && numnodes != 10) ExceptionT::GeneralFail(caller);
	if (numint != 1 &&
	    numint != 4) ExceptionT::GeneralFail(caller);	
	
	/* initialize */
	extrap = 0.0;

	switch (numint)
	{
		case 1:	
			
		extrap = 1.0;
		break;

		case 4:	
		{		
			double dat[4*10] = {
 1.9270509662496849, -0.3090169887498948, -0.3090169887498948, -0.3090169887498948, 
 0.809016988749895,  -0.3090169887498948,  0.809016988749895,   0.809016988749895, 
-0.3090169887498948, -0.3090169887498948, -0.3090169887498949,  1.9270509662496846,
-0.3090169887498949, -0.3090169887498949,  0.8090169887498949,  0.8090169887498949, 
-0.3090169887498949, -0.3090169887498949,  0.8090169887498949, -0.3090169887498949, 
-0.309016988749895,  -0.309016988749895,   1.9270509662496846, -0.309016988749895,
-0.309016988749895,   0.8090169887498948,  0.8090169887498948, -0.309016988749895, 
-0.309016988749895,   0.8090169887498948, -0.3090169887498949, -0.3090169887498949, 
-0.3090169887498949,  1.9270509662496846, -0.3090169887498949, -0.3090169887498949,
-0.3090169887498949,  0.8090169887498949,  0.8090169887498949,  0.8090169887498949
			};	

			dMatrixT smooth(10, 4, dat);
			smooth.CopyBlock(0, 0, extrap);
			extrap = smooth;
			
			break;
		}	
		default:
			ExceptionT::GeneralFail(caller);
	}
}

/* integration point gradient matrix */
void TetrahedronT::IPGradientTransform(int ip, dMatrixT& transform) const
{
	const char caller[] = "TetrahedronT::IPGradientTransform";

	/* dimensions */
	int nsd = transform.Rows();
	int nip = transform.Cols();
	if (nsd != 3) ExceptionT::SizeMismatch(caller);
	
	//TEMP only implemented for 1 integration point
	if (ip != 0 || nip != 1) ExceptionT::GeneralFail(caller, "only implemented for 1 integration point");

	/* no gradient */
	transform = 0.0;
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void TetrahedronT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	const char caller[] = "TetrahedronT::NodesOnFacet";

	if (fNumNodes != 4 && fNumNodes != 10)
		ExceptionT::GeneralFail(caller, "only implemented 4 and 10 element nodes: %d", fNumNodes);

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 4) ExceptionT::OutOfRange(caller);
#endif

	/* nodes-facet data */
	int dat4[] = {0,1,3,
		      1,2,3,
		      2,0,3,
		      0,2,1};
	
	int dat10[] = {0,1,3,4,8,7,
		       1,2,3,5,9,8,
		       2,0,3,6,7,9,
		       0,2,1,6,5,4};

	/* collect facet data */		
	iArrayT tmp;
	if (fNumNodes == 4)
		tmp.Set(3, dat4 + facet*3);
	else
		tmp.Set(6, dat10 + facet*6);
	
	/* (allocate and) copy in */
	facetnodes = tmp;
}

void TetrahedronT::NumNodesOnFacets(iArrayT& num_nodes) const
{
	if (fNumNodes != 4 && fNumNodes != 10)
		ExceptionT::GeneralFail("TetrahedronT::NodesOnFacet", "only implemented 4 and 10 element nodes");

	num_nodes.Dimension(4);
	if (fNumNodes == 4)
		num_nodes = 3;
	else
		num_nodes = 6;
}

/* return the local node numbers for each edge of element */
void TetrahedronT::NodesOnEdges(iArray2DT& nodes_on_edges) const
{
	/* nodes in edges data */
	int dat3[6*2] = {
		0,1,
		1,2,
		2,0,
		0,3,
		1,3,
		2,3
	};

	int dat10[6*3] = {
		0,4,1,
		1,5,2,
		2,6,0,
		0,7,3,
		1,8,3,
		2,10,3
	};

	iArray2DT tmp;
	if (fNumNodes == 3)
		tmp.Alias(6, 2, dat3);
	else if (fNumNodes == 10)
		tmp.Alias(6, 3, dat10);
	else
		ExceptionT::OutOfRange("TetrahedronT::NodesOnEdges");

	/* copy in */
	nodes_on_edges = tmp;
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void TetrahedronT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	int dat4[] = {0,1,3,
	              1,2,3,
	              2,0,3,
	              0,2,1};
	iArray2DT temp(4, 3, dat4);

	facetnodes = temp;
}

/* return geometry and number of nodes on each facet */
void TetrahedronT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
{
	if (fNumNodes != 4 && fNumNodes != 10)
		ExceptionT::GeneralFail("TetrahedronT::FacetGeometry", "only implemented for 4 nodes: %d", fNumNodes);

	facet_geom.Dimension(fNumFacets);
	facet_geom = kTriangle;
	
	facet_nodes.Dimension(fNumFacets);
	if (fNumNodes == 4)
		facet_nodes = 3;
	else
		facet_nodes = 6;
}

/* return true if the given point is within the domain */
bool TetrahedronT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
#if __option(extended_errorcheck)
		if (coords.NumberOfNodes() != 4) 
			ExceptionT::GeneralFail("TetrahedronT::PointInDomain", "expecting 4 element nodes: %d", coords.NumberOfNodes());
#endif

	/* nodes-facet data - ordered for outward normals */
	int dat4[] = {0,1,3,
	              1,2,3,
	              2,0,3,
	              0,2,1};

	/* method: check all faces and see of point lies inside */
	bool in_domain = true;
	int* facet_nodes = dat4;
	for (int i = 0; in_domain && i < 4; i++)
	{
		double ab_0 = coords(facet_nodes[1], 0) - coords(facet_nodes[0], 0);
		double ab_1 = coords(facet_nodes[1], 1) - coords(facet_nodes[0], 1);
		double ab_2 = coords(facet_nodes[1], 2) - coords(facet_nodes[0], 2);
		double ab_max = (fabs(ab_0) > fabs(ab_1)) ? fabs(ab_0) : fabs(ab_1);
		ab_max = (fabs(ab_2) > ab_max) ? fabs(ab_2) : ab_max;

		double ac_0 = coords(facet_nodes[2], 0) - coords(facet_nodes[0], 0);
		double ac_1 = coords(facet_nodes[2], 1) - coords(facet_nodes[0], 1);
		double ac_2 = coords(facet_nodes[2], 2) - coords(facet_nodes[0], 2);
		double ac_max = (fabs(ac_0) > fabs(ac_1)) ? fabs(ac_0) : fabs(ac_1);
		ac_max = (fabs(ac_2) > ac_max) ? fabs(ac_2) : ac_max;

		double L_ref = (ab_max > ac_max) ? ab_max : ac_max;

		double ap_0 = point[0] - coords(facet_nodes[0], 0);
		double ap_1 = point[1] - coords(facet_nodes[0], 1);
		double ap_2 = point[2] - coords(facet_nodes[0], 2);
			
		/* vector triple product */
		double ac_ab_0 = ac_1*ab_2 - ac_2*ab_1;
		double ac_ab_1 = ac_2*ab_0 - ac_0*ab_2;
		double ac_ab_2 = ac_0*ab_1 - ac_1*ab_0;			
		double triple_product = ac_ab_0*ap_0 + ac_ab_1*ap_1 + ac_ab_2*ap_2;
		in_domain = (triple_product/(L_ref*L_ref*L_ref)) > -kSmall;

		/* next face */		
		facet_nodes += 3;
	}
	
	return in_domain;
}
