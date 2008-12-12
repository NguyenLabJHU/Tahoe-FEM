/* $Id: GeometryBaseT.h,v 1.14 2008-12-12 17:44:27 lxmota Exp $ */
/* created: paklein (10/21/1997) */
#ifndef _GEOMETRY_BASE_T_H_
#define _GEOMETRY_BASE_T_H_

/* base class */
#include "GeometryT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;
class iArrayT;
class dArrayT;
class dMatrixT;
class iArray2DT;
template <class TYPE> class ArrayT;
class LocalArrayT;

/** base class for parent domain geometries. Derived classes must
 * initialize shape function arrays with geometry specific values. */
class GeometryBaseT: public GeometryT
{
public:

	/** constructor */
	GeometryBaseT(int numnodes, int numfacets);

	/** destructor */
	virtual ~GeometryBaseT(void);

	virtual const dArray2DT& ParentCoords(void) const;
  virtual const LocalArrayT& ParamCoords(void) const;
  virtual void SetParamCoords();

	/** returns the number of element facets */
	int NumFacets(void) const { return fNumFacets; };

	/** returns the number of element nodes */
	int NumNodes(void) const { return fNumNodes; };

	/** return the geometry code */
	virtual GeometryT::CodeT Geometry(void) const = 0;

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

	/** evaluate the shape functions and their first and second derivatives. the second
	 * derivative has been implemented for 27 node element only. */
	virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na,
		dArray2DT& DNa, dArray2DT& DDNa) const = 0;

	/** compute local shape functions and derivatives. The shape functions
	 * and their derivatives are evaluated for one of the pre-defined
	 * integration rules. The integration rule is determined from the
	 * dimensions of the arrays passed in.
	 * \param Na destination for shape function evaluations: [nip] x [nnd]
	 * \param Na_x destination for shape function gradients: [nip] x [nsd] x [nnd]
	 * \param weights destination for weights of the integration rule: [nip] */
	virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
		dArrayT& weights) const = 0;

	/* this function will be called from ParenDomainT.cpp to initialize
	 * local second derivative of shape functions.
	 * Local derivative of shape functions has been implemented for 27 node hex only */
	virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
		ArrayT<dArray2DT>& Na_xx, dArrayT& weights) const = 0;

	/** compute gradients of the "bubble" modes */
	virtual void BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const;

	/** set the values of the nodal extrapolation matrix */
	virtual void SetExtrapolation(dMatrixT& extrap) const = 0;

	/** integration point gradient matrix. Returns the matrix which transforms
	 * a vector of integration point values of a field to the gradient of that
	 * field wrt to the parent domain coordinates at the specified integration
	 * point. The matrix is define such that the gradient is given by
	 \f[
	 	\frac{\partial f}{\partial \xi_i} = A_{iI} f_I
	 \f]
	 * where \f$ f_I \f$ is a vector of the integration point values.
	 * This method must be overridden by derived types. GeometryBaseT::IPGradientTransform
	 * throws an exception if called. */
	virtual void IPGradientTransform(int ip, dMatrixT& transform) const;

	/** return the local node numbers for each facet of the element
	 * numbered to produce at outward normal in the order: vertex
	 * nodes, mid-edge nodes, mid-face nodes */
	virtual void NodesOnFacet(int facet, iArrayT& facets) const = 0;
	virtual void NumNodesOnFacets(iArrayT& num_nodes) const = 0;
//	virtual void NodesOnFacet(RaggedArrayT<int>& facets) const = 0;

	/** return the local node numbers for each edge of element where
	 * the row number corresponds with the canonical numbering of the
	 * edges in the element. Edges are only defined for 3D domain geometries. */
	virtual void NodesOnEdges(iArray2DT& nodes_on_edges) const = 0;

	/** returns the nodes on each facet needed to determine neighbors
	 * across facets. The number of nodes on each face may be less
	 * than the total number of nodes on a face if those additional
	 * nodes are redundant for determining the element neighbors, i.e.,
	 * mid-side nodes. */
	virtual void NeighborNodeMap(iArray2DT& facetnodes) const = 0;

	/** return geometry and number of nodes on each facet */
	virtual void FacetGeometry(ArrayT<CodeT>& facet_geom,
		iArrayT& facet_nodes) const = 0;

	/** return true if the given point is within the domain defined by
	 * the list of coordinates. Method must be overriden by subclasses.
	 * GeometryBaseT::PointInDomain throws an exception.
	 * \param coords list of coordinates defining the domain
	 * \param point test point coordinates */
	virtual bool PointInDomain(const LocalArrayT& coords, const dArrayT& point) const;

	/** return the integration point whose domain contains the given point in the
	 * parent domain coordinates */
	virtual int IPDomain(int nip, const dArrayT& coords) const;

	/** \name nodal subdomains
	 * Methods for accessing information about the nodal subdomains. These are the
	 * domains within the parent domain associated with each of the nodes. The
	 * union of subdomains from all parent domains comprising the support of a given node
	 * define a closed volume around the node that does not intersect the domain of
	 * any other node. */
	/*@{*/
	/** subdomain geometry */
	virtual GeometryT::CodeT NodalSubDomainGeometry(void) const;

	/** number of nodes defining the nodal subdomain */
	virtual int NodalSubDomainNumPoints(void) const;

	/** compute the coordinates of the points defining the nodal subdomain
	 * \param coords coordinates of the nodes over the entire parent domain
	 * \param node \e local number of node for which to compute the subdomain
	 * \param subdomain_coords returns with the coordinates of the points
	 *        defining the subdomain of the given node. */
	virtual void NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
		LocalArrayT& subdomain_coords) const;
	/*@}*/

protected:

	/* number of domain nodes */
	int fNumNodes;
	int fNumFacets;

	/*canonical coordinates of element nodes*/
	dArray2DT fCoords;
	LocalArrayT fParamCoords;
};

inline const dArray2DT& GeometryBaseT::ParentCoords(void) const
{
	return(fCoords);
}

inline const LocalArrayT& GeometryBaseT::ParamCoords() const
{
  return fParamCoords;
}

inline void GeometryBaseT::SetParamCoords() {

  fParamCoords.Dimension(fCoords.MinorDim(), fCoords.MajorDim());
  fParamCoords.Copy(fCoords.MinorDim(), fCoords.MajorDim(), fCoords);

}

} // namespace Tahoe
#endif /* _GEOMETRY_BASE_T_H_ */
