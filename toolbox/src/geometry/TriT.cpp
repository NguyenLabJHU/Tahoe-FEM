/* $Id: TriT.cpp,v 1.8 2005-03-02 02:27:14 paklein Exp $ */
/* created: paklein (07/03/1996) */
#include "TriT.h"
#include "QuadT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dMatrixT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* parameters */
const int kTrinsd           = 2;
const int kNumVertexNodes	= 3;

/* constructor */
TriT::TriT(int numnodes): GeometryBaseT(numnodes, kNumVertexNodes) {}

/* evaluate the shape functions and gradients. */
void TriT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
	const char caller[] = "TriT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
	if (coords.Length() != kTrinsd ||
	        Na.Length() != fNumNodes) ExceptionT::SizeMismatch(caller);
	if (fNumNodes != kNumVertexNodes) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
#if 0
	if (r < 0.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < 0.0 || s > 1.0) ExceptionT::OutOfRange(caller);
#endif

	/* shape functions */
	Na[0] = r;
	Na[1] = s;
	Na[2] = 1 - r - s;
}

/* evaluate the shape functions and gradients. */
void TriT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
{
	const char caller[] = "TriT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
	if (coords.Length() != kTrinsd ||
	        Na.Length() != fNumNodes ||
	     DNa.MajorDim() != kTrinsd ||
	     DNa.MinorDim() != fNumNodes) ExceptionT::SizeMismatch(caller);
	if (fNumNodes != kNumVertexNodes) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
#if 0
	if (r < 0.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < 0.0 || s > 1.0) ExceptionT::OutOfRange(caller);
#endif

	/* shape functions */
	Na[0] = r;
	Na[1] = s;
	Na[2] = 1 - r - s;

	/* derivatives */
	double* nax = DNa(0);
	double* nay = DNa(1);

	/* Na,r */
	nax[0] = 1.0;
	nax[1] = 0.0;
	nax[2] =-1.0;
	
	/* Na,s */
	nay[0] = 0.0;
	nay[1] = 1.0;
	nay[2] =-1.0;
}

/* compute local shape functions and derivatives */
void TriT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
	const char caller[] = "TriT::SetLocalShape";

	/* dimensions */
	int numnodes  = Na.MinorDim();
	int numint    = weights.Length();
	int nsd       = Na_x[0].MajorDim();

	/* dimension checks */
	if (numnodes < 3 || numnodes > 6)
		ExceptionT::GeneralFail(caller, "unsupported number of element nodes: %d", numnodes);

	if (numint != 1 &&
	    numint != 4 &&
	    numint != 6)
		ExceptionT::GeneralFail(caller, "unsupported number of integration points: %d", numint);

	if (nsd != kTrinsd) ExceptionT::GeneralFail(caller);

	/* initialize */
	Na = 0.0;
	for (int j = 0; j < Na_x.Length(); j++)
		Na_x[j] = 0.0;

	/* integration point coordinates */
	
	/* 1 point */
	double r1[1] = {1.0/3.0};
	double s1[1] = {1.0/3.0};
	
	/* 4 point */
	double r4[4] = {0.6, 0.2, 0.2, 1.0/3.0};
	double s4[4] = {0.2, 0.6, 0.2, 1.0/3.0};
	
	/* 6 point */
	double a1 = 0.816847572980459;
	double a2 = 0.091576213509771;
	double b1 = 0.108103018168070;
	double b2 = 0.445948490915965;
	double r6[6] = {a1, a2, a2, b1, b2, b2};
	double s6[6] = {a2, a1, a2, b2, b1, b2};
	
	double* r;
	double* s;

	/* set weights and r,s */
	switch (numint)
	{
		case 1:	
			
		weights[0] = 0.5; /* single integration point at the centroid */

		/* set coordinates */
		r = r1;
		s = s1;
		
		break;

		case 4:

		weights[0] = 25.0/48.0;
		weights[1] = weights[0];
		weights[2] = weights[0];
		weights[3] =-0.56250;  //centroid point
		weights *= 0.5;
		
		/* set coordinates */
		r = r4;
		s = s4;
		
		break;

	case 6:
	
		weights[0] = weights[1] = weights[2] = 0.109951743655322;
		weights[3] = weights[4] = weights[5] = 0.223381589678011;
		weights *= 0.5;
		
		/* set coordinates */
		r = r6;
		s = s6;
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

		/* vertex nodes */

	/* Na */
	na[0] += r[i];
	na[1] += s[i];
	na[2] += 1 - r[i] - s[i];

	/* Na,r */
	nax[0] += 1.0;
	nax[1] += 0.0;
	nax[2] +=-1.0;
	
	/* Na,s */
	nay[0] += 0.0;
	nay[1] += 1.0;
	nay[2] +=-1.0;
	
	/* mid-side nodes */
	if (numnodes > kNumVertexNodes)
	{
		/* add node 4 */
		if (numnodes > 3)
		{
			na[3] = 4.0*r[i]*s[i];	

			nax[3] = 4.0*s[i];
			nay[3] = 4.0*r[i];
		}
			
		/* add node 5 */
		if (numnodes > 4)
		{
			na[4] = 4.0*(1.0 - r[i] - s[i])*s[i];

			nax[4] =-4.0*s[i];
			nay[4] = 4.0*(1.0 - r[i] - 2.0*s[i]);
		}
			
		/* add node 6 */
		if (numnodes > 5)
		{
			na[5] = 4.0*(1.0 - r[i] - s[i])*r[i];			

			nax[5] = 4.0*(1.0 - 2.0*r[i] - s[i]);
			nay[5] =-4.0*r[i];
		}

			int corners[3][2] = {{0,1},	//adjacent to node 4
			                     {1,2}, //adjacent to node 5
			                     {2,0}};//adjacent to node 6

			/* corrections to vertex nodes */
			int mid = 0;
			for (int j = kNumVertexNodes; j < numnodes; j++) //mid-sides
			{
				for (int k = 0; k < 2; k++) //adjacent corners				
				{
	    			na[corners[mid][k]] -= 0.5*na[j];

	     			nax[corners[mid][k]] -= 0.5*nax[j];
	    			nay[corners[mid][k]] -= 0.5*nay[j];
				}
				
				mid++;
			}
	}	
	}
}

/* set the values of the nodal extrapolation matrix */
void TriT::SetExtrapolation(dMatrixT& extrap) const
{
	/* dimensions */
	int numnodes = extrap.Rows();
	int numint   = extrap.Cols();

	/* dimension checks */
	if (numnodes < 3 || numnodes > 6) throw ExceptionT::kGeneralFail;
	if (numint != 1 &&
	    numint != 4 &&
	    numint != 6) throw ExceptionT::kGeneralFail;	
	
	/* initialize */
	extrap = 0.0;

	switch (numint)
	{
		case 1:	
			
		extrap = 1.0;
		break;

		case 4:	

			/* vertex nodes */			
		extrap(0,0) = 5.0;
		extrap(0,1) = 2.5;
		extrap(0,2) = 2.5;
		extrap(0,3) =-9.0;

		extrap(1,0) = 2.5;
		extrap(1,1) = 5.0;
		extrap(1,2) = 2.5;
		extrap(1,3) =-9.0;

		extrap(2,0) =-1.25;
		extrap(2,1) =-1.25;
		extrap(2,2) = 1.25;
		extrap(2,3) = 2.25;

			/* mid-side nodes */
			if (numnodes > kNumVertexNodes)
			{
				double smooth[3][4] = {{-0.9375,-0.9375,-2.1875, 5.0625},
				                       { 0.625 , 1.875 , 1.875 ,-3.375 },
				                       { 1.875 , 0.625 , 1.875 ,-3.375 }};
			
				for (int i = kNumVertexNodes; i < numnodes; i++)
					for (int j = 0; j < numint; j++)
						extrap(i,j) = smooth[i-kNumVertexNodes][j];
			}
			
			break;

	case 6:
	{
		//Note: This nodal extrapolation rule was derived by fitting a
		//      bilinear surface to the values and locations of the
		//      centroids of each of the 4 sub-domains formed by the
		//      6 integration points.
	
		double smooth[6][6] =
		{{ 2.083114522569561, 1.142484133044889, 1.142484133044894,
		  -1.749781189236227,-0.809150799711556,-0.809150799711561},
{ 1.142484133044889, 2.083114522569561, 1.142484133044894,
-0.809150799711556,-1.749781189236226,-0.809150799711561},
{-0.5712420665224448,-0.5712420665224448, 0.3693883230022253,
0.904575399855778, 0.904575399855778,-0.03605498966889208},
{-0.3779681140117783,-0.3779681140117785,-0.848283308774116,
0.7113014473451123, 0.7113014473451123, 1.181616642107449},
{ 0.2856210332612223, 0.7559362280235579, 0.7559362280235597,
0.04771230007211126, -0.4226028946902242, -0.4226028946902264},
{0.7559362280235579, 0.2856210332612222, 0.7559362280235595,
-0.4226028946902244, 0.04771230007211102, -0.4226028946902267}};

//could really do them all like this.

			for (int i = 0; i < numnodes; i++)
				for (int j = 0; j < numint; j++)
					extrap(i,j) = smooth[i][j];

		break;
	}
	
		default:
		
			throw ExceptionT::kGeneralFail;
	}
}

/* integration point gradient matrix */
void TriT::IPGradientTransform(int ip, dMatrixT& transform) const
{
	const char caller[] = "TriT::IPGradientTransform";

	/* dimensions */
	int nsd = transform.Rows();
	int nip = transform.Cols();
	if (nsd != 2) ExceptionT::SizeMismatch(caller);
	
	//TEMP only implemented for 1 integration point
	if (ip != 0 || nip != 1) ExceptionT::GeneralFail(caller, "only implemented for 1 integration point");

	/* no gradient */
	transform = 0.0;
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void TriT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	const char caller[] = "TriT::NodesOnFacet";

	// TEMP: not implemented with midside nodes
	if (fNumNodes != 3 && fNumNodes != 6)
		ExceptionT::GeneralFail(caller, "only implemented for 3 and 6 element nodes: %d", fNumNodes);

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 2) ExceptionT::OutOfRange(caller);
#endif

	/* nodes-facet data */
	int dat3[] = {0,1,1,2,2,0};
	int dat6[] = {0,1,3,1,2,4,2,0,5};

	/* collect facet data */		
	iArrayT tmp;
	if (fNumNodes == 3)
		tmp.Set(2, dat3 + facet*2);
	else
		tmp.Set(3, dat6 + facet*3);
	
	/* (allocate and) copy in */
	facetnodes = tmp;
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void TriT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	int dat3[] = {0,1,1,2,2,0};
	iArray2DT temp(3, 2, dat3);
	
	facetnodes = temp;
}

void TriT::NumNodesOnFacets(iArrayT& num_nodes) const
{
// TEMP: not implemented with midside nodes
	if (fNumNodes != 3 && fNumNodes != 6)
		ExceptionT::GeneralFail("TriT::NumNodesOnFacets", "only implemented for 3 and 6 element nodes: %d", fNumNodes);

	num_nodes.Dimension(3);
	if (fNumNodes == 3)
		num_nodes = 2;
	else
		num_nodes = 3;
}

/* return geometry and number of nodes on each facet */
void TriT::FacetGeometry(ArrayT<CodeT>& facet_geom,
		iArrayT& facet_nodes) const
{
	facet_geom.Dimension(fNumFacets);
	facet_geom = kLine;
	
	facet_nodes.Dimension(fNumFacets);
	facet_nodes = 2;
	for (int i = 0; i < (fNumNodes - kNumVertexNodes); i++)
		facet_nodes[i] = 3;
}

/* return true if the given point is within the domain */
bool TriT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
	/* method: run around the perimeter of the element and see if
	 *         the point always lies to the left of segment a-b */
	int nen = coords.NumberOfNodes();
	int a = nen - 1;
	int b = 0;
	bool in_domain = true;
	for (int i = 0; in_domain && i < nen; i++)
	{
		double ab_0 = coords(b,0) - coords(a,0);
		double ab_1 = coords(b,1) - coords(a,1);

		double ap_0 = point[0] - coords(a,0);
		double ap_1 = point[1] - coords(a,1);
		
		double cross = ab_0*ap_1 - ab_1*ap_0;
		in_domain = cross >= 0.0;
		a++; 
		b++;
		if (a == nen) a = 0;
	}
	return in_domain;
}

/* subdomain geometry */
GeometryT::CodeT TriT::NodalSubDomainGeometry(void) const
{
	/* limited support */
	if (fNumNodes != 3)
		ExceptionT::GeneralFail("TriT::NodalSubDomainGeometry",
			"unsupported number of nodes %d", fNumNodes);

	return GeometryT::kQuadrilateral;
}

/* number of nodes defining the nodal subdomain */
int TriT::NodalSubDomainNumPoints(void) const
{
	/* limited support */
	if (fNumNodes != 3)
		ExceptionT::GeneralFail("TriT::NodalSubDomainGeometry",
			"unsupported number of nodes %d", fNumNodes);

	return 4;
}
	
/* compute the coordinates of the points defining the nodal subdomain */
void TriT::NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
	LocalArrayT& subdomain_coords) const
{
	const char caller[] = "TriT::NodalSubDomainCoordinates";

#if __option(extended_errorcheck)
	/* limited support */
	if (fNumNodes != 3)
		ExceptionT::GeneralFail(caller, "unsupported number of nodes %d", fNumNodes);

	/* checks */
	if (coords.NumberOfNodes() != fNumNodes || 
		coords.MinorDim() != 2 ||
		node < 0 || node >= fNumNodes ||
		subdomain_coords.MinorDim() != 2 ||
		subdomain_coords.NumberOfNodes() != TriT::NodalSubDomainNumPoints())
		ExceptionT::SizeMismatch(caller);
#endif

	/* tri domain */
	int next_node[3] = {1,2,0};
	int back_node[3] = {2,0,1};

	/* quad domain */
	int next_node_quad[4] = {1,2,3,0};
	int back_node_quad[4] = {3,0,1,2};
	int diag_node_quad[4] = {2,3,0,1};

	const double* px = coords(0);
	const double* py = coords(1);
	
	int back = back_node[node];
	int back_quad = back_node_quad[node];
	subdomain_coords(back_quad,0) = 0.5*(px[node] + px[back]);
	subdomain_coords(back_quad,1) = 0.5*(py[node] + py[back]);

	subdomain_coords(node,0) = px[node];
	subdomain_coords(node,1) = py[node];
	
	int next = next_node[node];
	int next_quad = next_node_quad[node];
	subdomain_coords(next_quad,0) = 0.5*(px[node] + px[next]);
	subdomain_coords(next_quad,1) = 0.5*(py[node] + py[next]);

	int diag = diag_node_quad[node];
	double one_third = 1.0/3.0;
	subdomain_coords(diag,0) = one_third*(px[0] + px[1] + px[2]);
	subdomain_coords(diag,1) = one_third*(py[0] + py[1] + py[2]);
}
