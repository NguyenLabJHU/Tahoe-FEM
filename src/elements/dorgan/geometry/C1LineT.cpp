/* $Id: C1LineT.cpp,v 1.2 2003-10-08 21:04:54 rdorgan Exp $ */
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
C1LineT::C1LineT(int numnodes): C1GeometryBaseT(numnodes, kNumVertexNodes) { }

/* evaluate the shape functions and gradients. */
void C1LineT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
        /* set shape functions */
        Na[0] = 0.5*(1.0 - coords[0]);
        Na[1] = 0.5*(1.0 + coords[0]);
}

/* evaluate the shape functions and gradients. */
void C1LineT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa) const
{
#pragma unused(coords)
#pragma unused(Na)
#pragma unused(DNa)
        cout << "\n C1LineT::EvaluateShapeFunctions: not implemented" << endl;
        throw ExceptionT::kGeneralFail;
}

/* compute local shape functions and derivatives */
void C1LineT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
        dArrayT& weights) const
{
        /* dimensions */
        int nnd = Na.MinorDim();
        int nip = weights.Length();
        int nsd = Na_x[0].MajorDim();

        /* dimension checks */
        if (nsd != kC1Linensd) throw ExceptionT::kGeneralFail;

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
                {
                        cout << "\n C1LineT::SetExtrapolation: unsupported number of integration points: " << nip << endl;
                        throw ExceptionT::kGeneralFail;
                }
        }
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void C1LineT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
#if __option(extended_errorcheck)
        if (facet != 0 && facet != 1) throw ExceptionT::kOutOfRange;
#else
#pragma unused (facet)        
#endif

        facetnodes.Dimension(1);
        facetnodes[0] = facet;
}

void C1LineT::NumNodesOnFacets(iArrayT& num_nodes) const
{
        num_nodes.Dimension(2);
        num_nodes = 1;
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void C1LineT::NeighborNodeMap(iArray2DT& facetnodes) const
{
        facetnodes.Dimension(2,1);
        facetnodes(0,0) = 0;
        facetnodes(1,0) = 1;
}

/* return geometry and number of nodes on each facet */
void C1LineT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
{
        facet_geom.Dimension(fNumFacets);
        facet_geom = kPoint;
        
        facet_nodes.Dimension(fNumFacets);
        facet_nodes = 1;
}
