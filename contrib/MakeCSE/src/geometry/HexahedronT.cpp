/*
 * File: HexahedronT.cpp
 *
 */

/*
 * created      : PAK (10/22/97)
 * last modified: PAK (05/06/99)
 */

#include <math.h>

#include "ExceptionCodes.h"

#include "HexahedronT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

/* parameters */
const int kHexnsd         = 3;
const int kNumVertexNodes = 8;
const int kNumFacets      = 6;

const double sqrt3 = sqrt(3.0);

/* constructor */
HexahedronT::HexahedronT(int numnodes): GeometryBaseT(numnodes, kNumFacets) {}

/* compute local shape functions and derivatives */
void HexahedronT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights)
{
	/* dimensions */
	int numnodes  = Na.MinorDim();
	int numint    = weights.Length();
	int nsd       = Na_x[0].MajorDim();

	/* dimension checks */
	if (numnodes != 8) throw eGeneralFail;
	if (numint != 1 && 
	    numint != 8 &&
	    numint != 27 &&
	    numint != 64) throw eGeneralFail;
	if (nsd != kHexnsd) throw eGeneralFail;

	/* initialize */
	Na = 0.0;
	for (int i = 0; i < Na_x.Length(); i++)
		Na_x[i] = 0.0;

	/* integration point coordinates */
	double	ra[] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
	double  sa[] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
	double  ta[] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};
	double xa_64[64], ya_64[64], za_64[64];
	double *xa, *ya, *za;
	double 	g;
	
	/* integration weights */
	switch (numint)
	{
		case 1:	
			  
    		g = 0.0;
    		xa = ra;
    		ya = sa;
    		za = ta;
    		weights[0] = 8.0;
    		break;
    
  		case 8:
  	
  			g = 1.0/sqrt3;
    		xa = ra;
    		ya = sa;
    		za = ta;
       		weights = 1.0;
			break;

		case 27:
		{
			/* coordinates */
			double b1 = sqrt(3.0/5.0);
			double b_1D[3] = {-b1, 0.0, b1}; 
			
			/* weights */
			double w1 = 5.0/9.0;
			double w2 = 8.0/9.0;
			double w_1D[3] = {w1, w2, w1};
			int x_i = 0;
			int y_i = 0;
			int z_i = 0;
			for (int i = 0; i < 27; i++)
			{
				xa_64[i]   = b_1D[x_i];
				ya_64[i]   = b_1D[y_i];
				za_64[i]   = b_1D[z_i];
				weights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];
							
				if (++x_i == 3)
				{
					x_i = 0;
					if (++y_i == 3)
					{
						y_i = 0;
						z_i++;
					}
				}
			}						
		
			xa = xa_64;
			ya = ya_64;
			za = za_64;
			g  = 1.0;		
			break;
			
		}	
		case 64:
		{
			/* coordinates */
			double b1 = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
			double b2 = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
			double b_1D[4] = {-b2,-b1, b1, b2};

			/* weights */					
			double w1 = (18.0 + sqrt(30.0))/36.0; 
			double w2 = (18.0 - sqrt(30.0))/36.0;
			double w_1D[4] = {w2, w1, w1, w2};

			int x_i = 0;
			int y_i = 0;
			int z_i = 0;
			for (int i = 0; i < 64; i++)
			{
				xa_64[i]   = b_1D[x_i];
				ya_64[i]   = b_1D[y_i];
				za_64[i]   = b_1D[z_i];
				weights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];
							
				if (++x_i == 4)
				{
					x_i = 0;
					if (++y_i == 4)
					{
						y_i = 0;
						z_i++;
					}
				}
			}						
			
			xa = xa_64;
			ya = ya_64;
			za = za_64;
			g  = 1.0;		
			break;
		}	
		default:
		
			throw eGeneralFail;
	}
				
	/* corner nodes only */
	for (int in = 0; in < numint; in++)	
	{
		double* na  = Na(in);
		double* nax = Na_x[in](0);
		double* nay = Na_x[in](1);
		double* naz = Na_x[in](2);

      	double r = g*xa[in];
      	double s = g*ya[in];
      	double t = g*za[in];
      	
		for (int lnd = 0; lnd < kNumVertexNodes; lnd++)
		{
      		double tempr1 = 1.0 + ra[lnd]*r;
      		double temps1 = 1.0 + sa[lnd]*s;
      		double tempt1 = 1.0 + ta[lnd]*t;
      		
      		*na++  = 0.125*tempr1*temps1*tempt1;
      		*nax++ = 0.125*ra[lnd]*temps1*tempt1;
      		*nay++ = 0.125*tempr1*sa[lnd]*tempt1;
      		*naz++ = 0.125*tempr1*temps1*ta[lnd];
       	}
 	}
}

/* set the values of the nodal extrapolation matrix */
void HexahedronT::SetExtrapolation(dMatrixT& extrap)
{
	/* dimensions */
	int numnodes = extrap.Rows();
	int numint   = extrap.Cols();

	/* dimension checks */
	if (numnodes != 8) throw eGeneralFail;

	/* initialize */
	extrap = 0.0;
	
	switch (numint)
	{
		case 1:	
			  
    		extrap(0,0) = 1.0;
    		break;
    
  		case 8:
		{
			double data[64] = {
	(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
	(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.,
	(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
	(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.,
	(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.,
	(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
	(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.,
	(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
	(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.,
	(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
	(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.,
	(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
	(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
	(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.,
	(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
	(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.};
	
			dMatrixT smooth(8,8,data);
			extrap = smooth;			

			break;
		}
		default:
				
			cout << "\n QuadT::SetExtrapolation: no nodal extrapolation with Gauss rule: ";
			cout << numint << endl;			
	}
}

/* return the local node numbers for each facet of the element 
 * numbered to produce at outward normal in the order: vertex
 * nodes, mid-edge nodes, mid-face nodes */
void HexahedronT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	if (fNumNodes != 8 && fNumNodes != 20)
	{
		cout << "\n HexahedronT::NodesOnFacet: only implemented 8 and 20 element nodes" << endl;
		throw eGeneralFail;
	}

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 5) throw eOutOfRange;
#endif

	/* nodes-facet data */
	int dat8[] = {0,3,2,1,
		      4,5,6,7,
		      0,1,5,4,
		      1,2,6,5,
		      2,3,7,6,
		      3,0,4,7};

	int dat20[] = {0,3,2,1,11,10, 9, 8,
		       4,5,6,7,12,13,14,15,
		       0,1,5,4, 8,17,12,16,
		       1,2,6,5, 9,18,13,17,
		       2,3,7,6,10,19,14,18,
		       3,0,4,7,11,16,15,19};

	/* collect facet data */		
	iArrayT tmp;
	if (fNumNodes == 8)
		tmp.Set(4, dat8 + facet*4);
	else
		tmp.Set(8, dat20 + facet*8);
	
	/* (allocate and) copy in */
	facetnodes = tmp;
}

void HexahedronT::NumNodesOnFacets(iArrayT& num_nodes) const
{
//TEMP
	if (fNumNodes != 8 && fNumNodes != 20)
	{
		cout << "\n HexahedronT::NumNodesOnFacets: only implemented 8 and 20 element nodes" << endl;
		throw eGeneralFail;
	}

	num_nodes.Allocate(6);
	if (fNumNodes == 8)
		num_nodes = 4;
	else
		num_nodes = 8;
}

/* returns the nodes on each facet needed to determine neighbors
 * across facets */
void HexahedronT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	/* nodes-facet data */
	int dat8[] = {0,3,2,1,
		          4,5,6,7,
		          0,1,5,4,
		          1,2,6,5,
		          2,3,7,6,
		          3,0,4,7};
	iArray2DT temp(6, 4, dat8);

	facetnodes = temp;
}

/* return geometry and number of nodes on each facet */
void HexahedronT::FacetGeometry(ArrayT<GeometryCode>& facet_geom, iArrayT& facet_nodes) const
{
//TEMP: reminder
	if (fNumNodes != 8 && fNumNodes != 20)
	{
		cout << "\n HexahedronT::FacetGeometry: only implemented for 8 and 20 nodes" << endl;
		throw eGeneralFail;
	}

	facet_geom.Allocate(fNumFacets);
	facet_geom = kQuadrilateral;
	
	facet_nodes.Allocate(fNumFacets);
	facet_nodes = 4;
}
