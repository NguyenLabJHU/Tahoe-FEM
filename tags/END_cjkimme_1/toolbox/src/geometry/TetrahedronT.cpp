/* $Id: TetrahedronT.cpp,v 1.3.8.2 2003-11-10 21:14:08 cjkimme Exp $ */
/* created: paklein (10/22/1996) */

#include "TetrahedronT.h"
#include "QuadT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

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
#pragma unused(coords)
#pragma unused(Na)

	cout << "\n TetrahedronT::EvaluateShapeFunctions: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* evaluate the shape functions and gradients. */
void TetrahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
{
#pragma unused(coords)
#pragma unused(Na)
#pragma unused(DNa)

	cout << "\n TetrahedronT::EvaluateShapeFunctions: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* compute local shape functions and derivatives */
void TetrahedronT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
	/* dimensions */
	int numnodes = Na.MinorDim();
	int numint   = weights.Length();
	int nsd      = Na_x[0].MajorDim();

	/* dimension checks */
	if (numnodes != 4)
	{
		cout << "\n TetrahedronT::SetLocalShape: unsupported number of element nodes: "
		     << numnodes << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	if (numint != 1 &&
	    numint != 4)
	{
		cout << "\n TetrahedronT::SetLocalShape: unsupported number of integration points: "
		     << numint << endl;
		throw ExceptionT::kGeneralFail;
	}

	if (nsd != kTetnsd) throw ExceptionT::kGeneralFail;

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
		
			throw ExceptionT::kGeneralFail;			
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
	/* dimensions */
	int numnodes = extrap.Rows();
	int numint   = extrap.Cols();

	/* dimension checks */
	if (numnodes != 4) throw ExceptionT::kGeneralFail;
	if (numint != 1 &&
	    numint != 4) throw ExceptionT::kGeneralFail;	
	
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
		
			throw ExceptionT::kGeneralFail;
	}
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void TetrahedronT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	if (fNumNodes != 4 && fNumNodes != 10)
	{
		cout << "\n TetrahedronT::NodesOnFacet: only implemented 4 and 10 element nodes" << endl;
		throw ExceptionT::kGeneralFail;
	}

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 4) throw ExceptionT::kOutOfRange;
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
	{
		cout << "\n TetrahedronT::NodesOnFacet: only implemented 4 and 10 element nodes" << endl;
		throw ExceptionT::kGeneralFail;
	}

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
	{
		cout << "\n TetrahedronT::FacetGeometry: only implemented for 4 nodes" << endl;
		throw ExceptionT::kGeneralFail;
	}

	facet_geom.Dimension(fNumFacets);
	facet_geom = kTriangle;
	
	facet_nodes.Dimension(fNumFacets);
	if (fNumNodes == 4)
		facet_nodes = 3;
	else
		facet_nodes = 6;
}
