/* $Id: LineT.cpp,v 1.10 2005-03-02 17:37:30 paklein Exp $ */
/* created: paklein (04/25/1999) */
#include "LineT.h"

#include <math.h>

#include "ExceptionT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dMatrixT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* parameters */
const int kLinensd = 1;
const int kNumVertexNodes = 2;
const double sqrt3 = sqrt(3.0);

/* constructor */
LineT::LineT(int numnodes): GeometryBaseT(numnodes, kNumVertexNodes) { }

/* evaluate the shape functions and gradients. */
void LineT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
	/* linear only */
	if (fNumNodes > 2)
		ExceptionT::GeneralFail("LineT::EvaluateShapeFunctions", "linear only");

	/* set shape functions */
	Na[0] = 0.5*(1.0 - coords[0]);
	Na[1] = 0.5*(1.0 + coords[0]);
}

/* evaluate the shape functions and gradients. */
void LineT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
{
	/* linear only */
	if (fNumNodes > 2)
		ExceptionT::GeneralFail("LineT::EvaluateShapeFunctions", "linear only");

	/* set shape functions */
	Na[0] = 0.5*(1.0 - coords[0]);
	Na[1] = 0.5*(1.0 + coords[0]);

	/* shape function derivatives */
	DNa(0,0) = -0.5*coords[0];
	DNa(0,1) = 0.5*coords[0];
}

/* compute local shape functions and derivatives */
void LineT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
	/* dimensions */
	int nnd = Na.MinorDim();
	int nip = weights.Length();
	int nsd = Na_x[0].MajorDim();

	/* dimension checks */
	if (nsd != kLinensd) ExceptionT::GeneralFail("LineT::SetLocalShape");

	/* initialize */
	Na = 0.0;
	for (int i = 0; i < nip; i++)
		Na_x[i] = 0.0;

	/* 1 point */
	double r1[1] = {0.0};
	
	/* 2 point */
	double x2    = sqrt(1.0/3.0);
double r2[2] = {-x2, x2};

	/* 3 point */
	double    x3 = sqrt(3.0/5.0);
double r3[3] = {-x3, 0.0, x3};

	/* 4 point */
	double x41 = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
	double x42 = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
	double r4[4] = {-x42,-x41, x41, x42};
	
	/* integration coordinates and weights */
	double* xa;
	switch (nip)
	{
		case 1:
		{
			xa = r1;			
		weights[0] = 2.0;
		break;
		}
		case 2:  	
		{
			xa = r2;			   			
		weights[0] = 1.0;
		weights[1] = 1.0;
			break;
		}
		case 3:
		{
			xa = r3;
			double a = 5.0/9.0;
			double b = 8.0/9.0;  			
		weights[0] = a;
		weights[1] = b;
		weights[2] = a;
			break;
		}
		case 4:
		{
			xa = r4;			    			

			double w1 = (18.0 + sqrt(30.0))/36.0;
			double w2 = (18.0 - sqrt(30.0))/36.0;
		weights[0] = w2;
		weights[1] = w1;
		weights[2] = w1;
		weights[3] = w2;
			break;
		}	
		default:		
			
			cout << "\n LineT::SetLocalShape: unsupported number of integration points: " << nip << endl;
			throw ExceptionT::kGeneralFail;
	}

	/* set shape functions and derivatives */
	switch (nnd)
	{
		case 2:  	
		{
			for (int i = 0; i < nip; i++)	
			{
				double* na  = Na(i);
				double* nax = Na_x[i](0);

		    	/* Na */
		    	na[0] = 0.5*(1.0 - xa[i]);
		    	na[1] = 0.5*(1.0 + xa[i]);
		
		        /* Na,x */
		    	nax[0] =-0.5;
		    	nax[1] = 0.5;
		    }
			break;
		}
		case 3:
		{
			for (int i = 0; i < nip; i++)	
			{
				double* na  = Na(i);
				double* nax = Na_x[i](0);

		    	/* Na */
		    	na[0] =-xa[i]*0.5*(1.0 - xa[i]);
		    	na[1] = xa[i]*0.5*(1.0 + xa[i]);
		    	na[2] = (1.0 - xa[i])*(1.0 + xa[i]);
		
		        /* Na,x */
		    	nax[0] =-0.5 + xa[i];
		    	nax[1] = 0.5 + xa[i];
		    	nax[2] =-2.0*xa[i];
		    }
			break;
		}
		default:
		
			cout << "\n LineT::SetLocalShape: unsupported number of nodes: " << nnd << endl;
			throw ExceptionT::kGeneralFail;
	}
}

/* set the values of the nodal extrapolation matrix */
void LineT::SetExtrapolation(dMatrixT& extrap) const
{
	/* dimensions */
	int nnd = extrap.Rows();
	int nip = extrap.Cols();

	/* dimension checks */
	if (nnd != 2 &&
	    nnd != 3) throw ExceptionT::kGeneralFail;

	/* initialize */
	extrap = 0.0;
	
	switch (nip)
	{
		case 1:	
			
		extrap = 1.0;
		break;

		case 2:	
		{
			double dat_2[3*2] = {
1.36602540378,
-0.366025403784,
0.5,
-0.366025403784,
1.36602540378,
0.5};
			dMatrixT extrap_2(3, 2, dat_2);
			extrap_2.CopyBlock(0, 0, extrap);
		break;
		}
		case 3:	
		{
			double dat_3[3*3] = {
1.4788305577,
0.187836108965,
0.0,
-0.666666666667,
-0.666666666667,
1.0,
0.187836108965,
1.4788305577,
0.0};
			dMatrixT extrap_3(3, 3, dat_3);
			extrap_3.CopyBlock(0, 0, extrap);
		break;
		}
		case 4:	
		{
			double dat_4[3*4] = {
1.52678812546,
-0.113917196282,
-0.0923265984407,
-0.813632449487,
0.400761520312,
0.592326598441,
0.400761520312,
-0.813632449487,
0.592326598441,
-0.113917196282,
1.52678812546,
-0.0923265984407};
			dMatrixT extrap_4(3, 4, dat_4);
			extrap_4.CopyBlock(0, 0, extrap);
		break;
		}
	default:
		ExceptionT::GeneralFail("LineT::SetExtrapolation", "unsupported integration rule %d", nip);
	}
}

/* integration point gradient matrix */
void LineT::IPGradientTransform(int ip, dMatrixT& transform) const
{
	const char caller[] = "LineT::IPGradientTransform";

	/* dimensions */
	int nsd = transform.Rows();
	int nip = transform.Cols();
	if (nsd != 1) ExceptionT::SizeMismatch(caller);

	/* constant gradient */
#pragma unused(ip)
	if (nip == 1)
		transform = 0.0;
	else if (nip == 2) {
		double a = sqrt(3.0)/2.0;
		transform(0,0) =-a;
		transform(0,1) = a;
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported number of integration points %d", nip);
}
/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void LineT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
#if __option(extended_errorcheck)
	if (facet != 0 && facet != 1) throw ExceptionT::kOutOfRange;
#else
#pragma unused (facet)	
#endif

	facetnodes.Dimension(1);
	facetnodes[0] = facet;
}

void LineT::NumNodesOnFacets(iArrayT& num_nodes) const
{
	num_nodes.Dimension(2);
	num_nodes = 1;
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void LineT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	facetnodes.Dimension(2,1);
	facetnodes(0,0) = 0;
	facetnodes(1,0) = 1;
}

/* return geometry and number of nodes on each facet */
void LineT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
{
	facet_geom.Dimension(fNumFacets);
	facet_geom = kPoint;
	
	facet_nodes.Dimension(fNumFacets);
	facet_nodes = 1;
}

/* return true if the given point is within the domain */
bool LineT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
#if __option(extended_errorcheck)
	if (coords.NumberOfNodes() != 2) 
		ExceptionT::GeneralFail("LineT::PointInDomain", "expecting only 2 points: %d", coords.NumberOfNodes());
#endif

	double v_01 = coords[1] - coords[0];
	double v_0p =  point[0] - coords[0];
	double v_1p =  point[0] - coords[1];
	return (v_1p/v_01) < kSmall && (-v_0p/v_01) < kSmall;
}

/* return the integration point whose domain contains the given point in the parent domain coordinates */
int LineT::IPDomain(int nip, const dArrayT& coords) const
{
	const char caller[] = "LineT::IPDomain";
	
	/* domain check */
	if (coords[0] < -1.0 || coords[0] > 1.0)
		ExceptionT::OutOfRange(caller, "{%g} outside domain", coords[0]);	
	
	if (nip == 1)
		return 0;
	else if (nip == 2) {
		if (coords[0] > 0.0)
			return 1;
		else
			return 0;
	}
	else
		ExceptionT::GeneralFail(caller, "%d integration points not supported", nip);

	/* dummy */
	return -1;
}
