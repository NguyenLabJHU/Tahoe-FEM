/* $Id: TetrahedronT.cpp,v 1.6 2004-05-12 22:20:15 paklein Exp $ */
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
	if (fNumNodes != kNumVertexNodes) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	if (r < 0.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < 0.0 || s > 1.0) ExceptionT::OutOfRange(caller);
	if (t < 0.0 || t > 1.0) ExceptionT::OutOfRange(caller);

	/* shape functions */
	Na[0] = r;
	Na[1] = s;
	Na[3] = t;
	Na[2] = 1 - r - s - t;
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
	if (fNumNodes != kNumVertexNodes) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	if (r < 0.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < 0.0 || s > 1.0) ExceptionT::OutOfRange(caller);
	if (t < 0.0 || t > 1.0) ExceptionT::OutOfRange(caller);

	/* shape functions */
	Na[0] = r;
	Na[1] = s;
	Na[3] = t;
	Na[2] = 1 - r - s - t;

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
	if (numnodes != 4)
		ExceptionT::GeneralFail(caller, "unsupported number of element nodes: %d", numnodes);
	
	if (numint != 1 &&
	    numint != 4)
		ExceptionT::GeneralFail(caller, "unsupported number of integration points: %d", numint);

	if (nsd != kTetnsd) ExceptionT::GeneralFail(caller);

	/* initialize */
	Na = 0.0;
	for (int j = 0; j < Na_x.Length(); j++)
		Na_x[j] = 0.0;

	/* integration point coordinates */
	
	/* 1 point */
	double r1[1] = {1.0/4.0};
	double s1[1] = {1.0/4.0};
	double t1[1] = {1.0/4.0};
	
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
			
		weights[0] = 1.0/6.0;

			/* set coordinates */
		r = r1;
		s = s1;
		t = t1;
		
		break;

		case 4:

		weights = (1.0/24.0);
		
			/* set coordinates */
		r = r4;
		s = s4;
		t = t4;
		
		break;

		default:
		
			ExceptionT::GeneralFail(caller);
	}	

	/* shape functions and derivatives */
	for (int i = 0; i < numint; i++)
	{
		double* na  = Na(i);
		double* nax = Na_x[i](0);
		double* nay = Na_x[i](1);
		double* naz = Na_x[i](2);

		/* vertex nodes */

	/* Na */
	na[0] += r[i];
	na[1] += s[i];
	na[3] += t[i];
	na[2] += 1 - r[i] - s[i] - t[i];

	/* Na,r */
	nax[0] += 1.0;
	nax[1] += 0.0;
	nax[3] += 0.0;
	nax[2] +=-1.0;
	
	/* Na,s */
	nay[0] += 0.0;
	nay[1] += 1.0;
	nay[3] += 0.0;
	nay[2] +=-1.0;

	/* Na,t */
	naz[0] += 0.0;
	naz[1] += 0.0;
	naz[3] += 1.0;
	naz[2] +=-1.0;    	
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
	if (numnodes != 4) ExceptionT::GeneralFail(caller);
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
			double dat[16] = {
	 1.92705096625,  -0.30901698875,  -0.30901698875,  -0.30901698875,
	-0.30901698875,   1.92705096625,  -0.30901698875,  -0.30901698875,
	-0.30901698875,  -0.30901698875,   1.92705096625,  -0.30901698875,
	-0.30901698875,  -0.30901698875,  -0.30901698875,   1.92705096625};
	
			dMatrixT smooth(4,4,dat);
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
		/* facet 1 */
		double ab_0 = coords(facet_nodes[1], 0) - coords(facet_nodes[0], 0);
		double ab_1 = coords(facet_nodes[1], 1) - coords(facet_nodes[0], 1);
		double ab_2 = coords(facet_nodes[1], 2) - coords(facet_nodes[0], 2);

		double ac_0 = coords(facet_nodes[2], 0) - coords(facet_nodes[0], 0);
		double ac_1 = coords(facet_nodes[2], 1) - coords(facet_nodes[0], 1);
		double ac_2 = coords(facet_nodes[2], 2) - coords(facet_nodes[0], 2);

		double ap_0 = point[0] - coords(facet_nodes[0], 0);
		double ap_1 = point[1] - coords(facet_nodes[0], 1);
		double ap_2 = point[2] - coords(facet_nodes[0], 2);
			
		/* vector triple product */
		double ac_ab_0 = ac_1*ab_2 - ac_2*ab_1;
		double ac_ab_1 = ac_2*ab_0 - ac_0*ab_2;
		double ac_ab_2 = ac_0*ab_1 - ac_1*ab_0;			
		double triple_product = ac_ab_0*ap_0 + ac_ab_1*ap_1 + ac_ab_2*ap_2;
		in_domain = triple_product >= 0.0;

		/* next face */		
		facet_nodes += 3;
	}
	
	return in_domain;
}
