/* $Id: C1LineT.h,v 1.3 2003-11-19 20:38:24 rdorgan Exp $ */
#ifndef _C1_LINE_T_H_
#define _C1_LINE_T_H_

/* base class */
#include "C1GeometryBaseT.h"

namespace Tahoe {

/** 1D line element */
class C1LineT: public C1GeometryBaseT
{
public:

        /** constructor */
        C1LineT(int numnodes);

        /** return the geometry code */
        virtual C1GeometryT::CodeT C1Geometry(void) const { return kC1Line; };

        /** evaluate the shape functions. See 
         * C1GeometryBaseT::EvaluateShapeFunctions for documentation */
        virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const;

        //        /** evaluate the shape functions and gradients. See 
        //         * C1GeometryBaseT::EvaluateShapeFunctions for documentation */
        //        virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, 
        //                dArrayT& DDNa) const;

        /** evaluate the shape functions and gradients. See 
         * C1GeometryBaseT::SetLocalShape for documentation */
        virtual void SetLocalShape(dArray2DT& Na, dArray2DT& Na_xx, dArrayT& weights) const;

        /* set the values of the nodal extrapolation matrix */
        virtual void SetExtrapolation(dMatrixT& extrap) const;

        //        /* return the local node numbers for each facet of the element
        //         * numbered to produce at outward normal in the order: vertex
        //         * nodes, mid-edge nodes, mid-face nodes */
        //        virtual void NodesOnFacet(int facet, iArrayT& facetnodes) const;
        ////        virtual void NodesOnFacet(RaggedArrayT<int>& facets) const;
        //        virtual void NumNodesOnFacets(iArrayT& num_nodes) const;

        //        /* returns the nodes on each facet needed to determine neighbors
        //         * across facets */
        //        virtual void NeighborNodeMap(iArray2DT& facetnodes) const;

        //        /* return geometry and number of nodes on each facet */
        //        virtual void FacetGeometry(ArrayT<CodeT>& facet_geom,
        //                iArrayT& facet_nodes) const;

        /** \name accessors */
        /*@{*/
        int NumNodes(void) const;
        /*@}*/

private:

	/* number of domain nodes */
	int fNumNodes;
};

/* inlines */

/* number of spatial dimensions */
inline int C1LineT::NumNodes(void) const { return fNumNodes; }

} // namespace Tahoe 
#endif /* _C1_LINE_T_H_ */
