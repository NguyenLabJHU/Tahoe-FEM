/* $Id: C1LineT.cpp,v 1.3 2003-11-19 20:38:24 rdorgan Exp $ */
#include "C1LineT.h"

#include <math.h>

#include "ExceptionT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dMatrixT.h"


using namespace Tahoe;

/* parameters */
const int kC1Linensd = 1;
const int kNumVertexNodes = 2;
const double sqrt3 = sqrt(3.0);

/* constructor */
C1LineT::C1LineT(int numnodes):
        C1GeometryBaseT(numnodes),  //, kNumVertexNodes),
	fNumNodes(numnodes)
{

}

/* evaluate the shape functions and gradients. */
void C1LineT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
        /* set shape functions */
        Na[0] =  0.25*(coords[0] - 1)*(coords[0] - 1)*(coords[0] + 2);
	Na[1] = -0.25*(coords[0] + 1)*(coords[0] + 1)*(coords[0] - 2);

	Na[2] =  0.25*(coords[0] - 1)*(coords[0] - 1)*(coords[0] + 1);
	Na[3] =  0.25*(coords[0] + 1)*(coords[0] + 1)*(coords[0] - 1);
}

        ///* evaluate the shape functions and gradients. */
        //void C1LineT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
        //{
        //#pragma unused(coords)
        //#pragma unused(Na)
        //#pragma unused(DNa)
        //        cout << "\n C1LineT::EvaluateShapeFunctions: not implemented" << endl;
        //        throw ExceptionT::kGeneralFail;
        //}

/* compute local shape functions and derivatives */
void C1LineT::SetLocalShape(dArray2DT& Na, dArray2DT& Na_xx,
        dArrayT& weights) const
{
        //        /* dimensions */
        //        int nnd = Na.MinorDim();
        //        int nip = weights.Length();
        //        int nsd = Na_x[0].MajorDim();

	/* dimensions */
	int nip = weights.Length();
	int nnd = NumNodes();
	int nsd = kC1Linensd;

	double nsd_d = log(1.0 * Na.MinorDim() / nnd) / log(2.0);

        /* dimension checks */
        if (nsd_d != nsd) throw ExceptionT::kGeneralFail;

        /* initialize */
        Na = 0.0;
        Na_xx = 0.0;

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
                        
                        cout << "\n C1LineT::SetLocalShape: unsupported number of integration points: " << nip << endl;
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
                                double* naxx = Na_xx(i);

                                /* Na */
				na[0] =  0.25*(xa[i] - 1)*(xa[i] - 1)*(xa[i] + 2);
				na[1] = -0.25*(xa[i] + 1)*(xa[i] + 1)*(xa[i] - 2);

				na[2] =  0.25*(xa[i] - 1)*(xa[i] - 1)*(xa[i] + 1);
				na[3] =  0.25*(xa[i] + 1)*(xa[i] + 1)*(xa[i] - 1);
                
                                /* Na,xx */
				naxx[0] =  1.5*xa[i];
				naxx[1] = -1.5*xa[i];

				naxx[2] =  0.5*(3*xa[i] - 1);
				naxx[3] =  0.5*(3*xa[i] + 1);
                        }
                        break;
                }
                default:
                
                        cout << "\n C1LineT::SetLocalShape: unsupported number of nodes: " << nnd << endl;
                        throw ExceptionT::kGeneralFail;
        }
}

/* set the values of the nodal extrapolation matrix */
void C1LineT::SetExtrapolation(dMatrixT& extrap) const
{
        /* dimensions */
        int nnd = extrap.Rows();
        int nip = extrap.Cols();

        /* dimension checks */
        if (nnd != 2) throw ExceptionT::kGeneralFail;

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
                {
                        cout << "\n C1LineT::SetExtrapolation: unsupported number of integration points: " << nip << endl;
                        throw ExceptionT::kGeneralFail;
                }
        }
}

        ///* return the local node numbers for each facet of the element
        //* numbered to produce at outward normal in the order: vertex
        //* nodes, mid-edge nodes, mid-face nodes */
        //void C1LineT::NodesOnFacet(int facet, iArrayT& facetnodes) const
        //{
        //#if __option(extended_errorcheck)
        //        if (facet != 0 && facet != 1) throw ExceptionT::kOutOfRange;
        //#else
        //#pragma unused (facet)        
        //#endif
        //
        //        facetnodes.Dimension(1);
        //        facetnodes[0] = facet;
        //}

        //void C1LineT::NumNodesOnFacets(iArrayT& num_nodes) const
        //{
        //        num_nodes.Dimension(2);
        //        num_nodes = 1;
        //}

        ///* returns the nodes on each facet needed to determine neighbors
        //* across facets */
        //void C1LineT::NeighborNodeMap(iArray2DT& facetnodes) const
        //{
        //        facetnodes.Dimension(2,1);
        //        facetnodes(0,0) = 0;
        //        facetnodes(1,0) = 1;
        //}

        ///* return geometry and number of nodes on each facet */
        //void C1LineT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
        //{
        //        facet_geom.Dimension(fNumFacets);
        //        facet_geom = kPoint;
        //        
        //        facet_nodes.Dimension(fNumFacets);
        //        facet_nodes = 1;
        //}
