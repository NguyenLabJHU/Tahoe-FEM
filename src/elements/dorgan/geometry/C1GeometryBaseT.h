/* $Id: C1GeometryBaseT.h,v 1.2 2003-10-08 21:04:54 rdorgan Exp $ */
#ifndef _C1_GEOMETRY_BASE_T_H_
#define _C1_GEOMETRY_BASE_T_H_

/* base class */
#include "C1GeometryT.h"

namespace Tahoe {

/* forward declarations */
class dArray2DT;
template <class TYPE> class ArrayT;
class iArrayT;
class dArrayT;
class dMatrixT;
class iArray2DT;
template <class TYPE> class ArrayT;

/** base class for parent domain geometries. Derived classes must 
 * initialize shape function arrays with geometry specific values. */
class C1GeometryBaseT: public C1GeometryT
{
public:

        /** constructor */
        C1GeometryBaseT(int numnodes, int numfacets);

        /** destructor */
        virtual ~C1GeometryBaseT(void);

        /** returns the number of element facets */
        int NumFacets(void) const { return fNumFacets; };

        /** returns the number of element nodes */
        int NumNodes(void) const { return fNumNodes; };
        
        /** return the geometry code */
        virtual C1GeometryT::CodeT C1Geometry(void) const = 0;

        /** evaluate the shape functions. Compute the values of the
         * shape functions at an arbirary point in the
         * in the parent domain. Coordinates must fall within the domain.
         * \param coords point in the parent domain
         * \param Na destination for shape function values for each of the domain
         *        nodes. Must be dimensioned: [nnd] */
        virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const = 0;
        
        /** evaluate the shape functions and gradients. Compute the values of the
         * shape functions and their gradients at an arbirary point in the
         * in the parent domain. Coordinates must fall within the domain.
         * \param coords point in the parent domain
         * \param Na destination for shape function values for each of the domain
         *        nodes. Must be dimensioned: [nnd]
         * \param DNa destination for shape function derivatives. Must be 
         *        dimensioned: [nsd] x [nnd] */
        virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, 
                dArray2DT& DNa) const = 0;

        /** compute local shape functions and derivatives. The shape functions
         * and their derivatives are evaluated for one of the pre-defined
         * integration rules. The integration rule is determined from the
         * dimensions of the arrays passed in.
         * \param Na destination for shape function evaluations: [nip] x [nnd]
         * \param Na_x destination for shape function gradients: [nip] x [nsd] x [nnd]
         * \param weights destination for weights of the integration rule: [nip] */
        virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
                dArrayT& weights) const = 0;

        /** compute gradients of the "bubble" modes */
        virtual void BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const;

        /** set the values of the nodal extrapolation matrix */
        virtual void SetExtrapolation(dMatrixT& extrap) const = 0;
        
        /* return the local node numbers for each facet of the element
         * numbered to produce at outward normal in the order: vertex
         * nodes, mid-edge nodes, mid-face nodes */
        virtual void NodesOnFacet(int facet, iArrayT& facets) const = 0;
        virtual void NumNodesOnFacets(iArrayT& num_nodes) const = 0;
//        virtual void NodesOnFacet(RaggedArrayT<int>& facets) const = 0;

        /** returns the nodes on each facet needed to determine neighbors
         * across facets. The number of nodes on each face may be less
         * than the total number of nodes on a face if those additional
         * nodes are redundant for determining the element neighbors, i.e.,
         * mid-side nodes. */
        virtual void NeighborNodeMap(iArray2DT& facetnodes) const = 0;

        /** return geometry and number of nodes on each facet */
        virtual void FacetGeometry(ArrayT<CodeT>& facet_geom,
                iArrayT& facet_nodes) const = 0;
                
protected:

        /* number of domain nodes */
        int fNumNodes;
        int fNumFacets;
};

} // namespace Tahoe 
#endif /* _C1_GEOMETRY_BASE_T_H_ */
