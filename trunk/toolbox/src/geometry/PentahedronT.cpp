/* $Id: PentahedronT.cpp,v 1.7 2006-10-26 19:07:04 regueiro Exp $ */
/* created: sawimme (10/22/1999) */
#include "PentahedronT.h"

#include "ExceptionT.h"
#include "QuadT.h"
#include "TriT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* parameters */
const int kTetnsd           = 3;
const int kNumVertexNodes   = 6;
const int kNumFacets        = 5;

/* constructor */
PentahedronT::PentahedronT(int numnodes): GeometryBaseT(numnodes, kNumFacets) 
{
}

const dArray2DT& PentahedronT::ParentCoords(void) const
{
	cout << "\n PentahedronT::ParentCoords: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* evaluate the shape functions and gradients. */
void PentahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
#pragma unused(coords)
#pragma unused(Na)

	cout << "\n PentahedronT::EvaluateShapeFunctions: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* evaluate the shape functions and gradients. */
void PentahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
{
#pragma unused(coords)
#pragma unused(Na)
#pragma unused(DNa)

	cout << "\n PentahedronT::EvaluateShapeFunctions: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}


/* evaluate the shape functions and first and second gradients. */
void PentahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa, dArray2DT& DDNa) const
{
#pragma unused(coords)
#pragma unused(Na)
#pragma unused(DNa)
#pragma unused(DDNa)

	cout << "\n PentahedronT::EvaluateShapeFunctions: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}


/* compute local shape functions and derivatives */
void PentahedronT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
#pragma unused(Na)
#pragma unused(Na_x)
#pragma unused(weights)

cout << "\n PentahedronT::SetLocalShape not implemented " << endl;
throw ExceptionT::kGeneralFail;
}

/* set the values of the nodal extrapolation matrix */
void PentahedronT::SetExtrapolation(dMatrixT& extrap) const
{
#pragma unused(extrap)

cout << "\n PentahedronT::SetExtrapolation not implemented " << endl;
throw ExceptionT::kGeneralFail;
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void PentahedronT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	if (fNumNodes != 6 && fNumNodes != 15)
	{
		cout << "\n PentahedronT::NodesOnFacet: only implemented 6 and 15 element nodes" << endl;
		throw ExceptionT::kGeneralFail;
	}

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 5) throw ExceptionT::kOutOfRange;
#endif

	/* nodes-facet data */
	int dat6_0to1[] = {0,2,1,
			   3,4,5};
	int dat6_2to5[] = {0,3,4,1,
			   1,4,5,2,
			   2,5,3,0};
	
	int dat15_0to1[] = {0,2,1, 6,8,7,
			    3,4,5, 11,9,10};
	int dat15_2to5[] = {0,3,4,1, 10,9,13,6,
			    1,4,5,2, 13,10,14,7,
			    2,5,3,0, 14,11,10,8};

	/* collect facet data */		
	iArrayT tmp;
	if (fNumNodes == 6)
	  {
	    if (facet < 2)
	      tmp.Set (3, dat6_0to1 + facet*3);
	    else
	      tmp.Set (4, dat6_2to5 + (facet - 2)*4);
	  }
	else
	  {
	    if (facet < 2)
	      tmp.Set (6, dat15_0to1 + facet*6);
	    else
	      tmp.Set (8, dat15_2to5 + (facet - 2)*8);
	  }
	
	/* (allocate and) copy in */
	facetnodes = tmp;
}

void PentahedronT::NumNodesOnFacets(iArrayT& num_nodes) const
{
	if (fNumNodes != 6 && fNumNodes != 15)
	{
		cout << "\n PentahedronT::NodesOnFacet: only implemented 6 and 15 element nodes" << endl;
		throw ExceptionT::kGeneralFail;
	}

	num_nodes.Dimension(5);
	if (fNumNodes == 6)
	{
		num_nodes = 4;
		num_nodes[0] = 3;
		num_nodes[1] = 3;
	}
	else
	{
		num_nodes = 8;
		num_nodes[0] = 6;
		num_nodes[1] = 6;
	}
}

/* return the local node numbers for each edge of element */
void PentahedronT::NodesOnEdges(iArray2DT& nodes_on_edges) const
{
	/* nodes in edges data */
	int dat6[9*2] = {
		0,1,
		1,2,
		2,0,
		3,4,
		4,5,
		5,3,
		0,3,
		1,4,
		2,5
	};

	int dat15[9*3] = {
		0,6,1,
		1,7,2,
		2,8,0,
		3,9,4,
		4,10,5,
		5,11,3,
		0,12,3,
		1,13,4,
		2,14,5
	};

	iArray2DT tmp;
	if (fNumNodes == 6)
		tmp.Alias(9, 2, dat6);
	else if (fNumNodes == 15)
		tmp.Alias(9, 3, dat15);
	else
		ExceptionT::OutOfRange("PentahedronT::NodesOnEdges");

	/* copy in */
	nodes_on_edges = tmp;
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void PentahedronT::NeighborNodeMap(iArray2DT& facetnodes) const
{
#pragma unused(facetnodes)

cout << "\n PentahedronT::NeighborNodeMap not implemented " << endl;
// this would require facetnodes to be a ragged array
throw ExceptionT::kGeneralFail;
}

/* return geometry and number of nodes on each facet */
void PentahedronT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
{
	facet_geom.Dimension(fNumFacets);
	facet_geom = kQuadrilateral;
	facet_geom[0] = kTriangle;
	facet_geom[1] = kTriangle;
	
	NumNodesOnFacets (facet_nodes);
}
