/* $Id: ParentDomainT.h,v 1.2 2001-06-02 02:20:43 paklein Exp $ */
/* created: paklein (07/03/1996)                                          */
/* interface for a finite element parent domain. manages integration      */
/* information (points, weights, etc..) and mapping between the real      */
/* coordinates and the domain.                                            */

#ifndef _PARENT_DOMAIN_T_H_
#define _PARENT_DOMAIN_T_H_

/* direct members */
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dArray2DT.h"
#include "GeometryBaseT.h"

/* forward declarations */
class iArrayT;
class LocalArrayT;

/** class to handle calculations over a parent domain */
class ParentDomainT
{
  public:
	/** constructor.
	 * \param geometry_code geometry of the domain
	 * \param numIP number of integration points to evaluate in the domain
	 * \param numnodes number of nodes defining the domain geometry */
	ParentDomainT(GeometryT::CodeT geometry_code, int numIP, int numnodes);

	/** destructor */
	~ParentDomainT(void);

	/** set all local parameters. call immediately after constructor */
	void Initialize(void);

	/* accessors */
	int NumSD(void) const;
	int NumIP(void) const;
	int NumNodes(void) const;
	GeometryT::CodeT GeometryCode(void) const;

	/** reference to the entire shape function array.
	 * \return 2D array: [nip] x [nnd] */
	const dArray2DT& Na(void) const;
	
	/** pointer to the shape functions.
	 * \param IPnum integration point number
	 * \return pointer to array length numnodes */
	const double* Shape(int IPnum) const;

	/** pointer to the shape functions derivatives.
	 * \param IPnum integration point number
	 * \param dim derivative component
	 * \return pointer to array length numnodes */
	const double* DShape(int IPnum, int dim) const;

	/** pointer to the weights for all integration points */
	const double* Weight(void) const;

	/** interpolation of nodal values.
	 * \param nodal values at the nodes
	 * \param interp result of the interpolation
	 * \param IPnum integration point number */
	void Interpolate(const LocalArrayT& nodal, dArrayT& interp, int IPnum) const;
	
	/** interpolation of nodal values to all integration points.
	 * \param nodal values at the nodes
	 * \param interp interpolation to all ip's: [nip] x [nu] */
	void Interpolate(const LocalArrayT& nodal, dArray2DT& interp) const;

	/** compute the jacobian of the nodal values.
	 * uses externally provided shape function derivatives.
	 * \param nodal values at the nodes: [nnd] x [nu]
	 * \param DNa shape function derivatives: [ndim] x [nnd]
	 * \param jacobian resulting jacobian: [nu] x [ndim] */
	void Jacobian(const LocalArrayT& nodal, const dArray2DT& DNa, dMatrixT& jacobian) const;

	/** compute the jacobian of the nodal values with respect to domain coordinates.
	 * \param nodal values at the nodes: [nnd] x [nu]
	 * \param numIP integration point number
	 * \param jacobian resulting jacobian: [nu] x [nsd] */
	void DomainJacobian(const LocalArrayT& nodal, int numIP, dMatrixT& jacobian) const;

	/** norm of the surface mapping.
	 * \param jacobian surface jacobian: [nsd] x [nsd - 1]
	 * \return jacobian of the surface transformation */
	double SurfaceJacobian(const dMatrixT& jacobian) const;

	/** surface transformations.
	 * compute the coordinate transformation and the norm of the
	 * surface mapping for the given jacobian.
	 * \param jacobian surface jacobian: [nsd] x [nsd - 1]
	 * \param  Q transformation from global to local(') coordinates, i.e., 
	 * t'_i = Q_ik t_k, where t'_j (j = nsd) is the "normal" direction
	 * \return jacobian of the surface transformation */
	double SurfaceJacobian(const dMatrixT& jacobian, dMatrixT& Q) const;

	/** shape function derivatives.
	 * compute the derivatives of the shape functions with respect
	 * to the given coordinates by the chain rule for all integration
	 * points at once.
	 * \param coords nodal coordinates
	 * \param DNa shape function derivatives with respect to given 
	 *        coordinates: [nip] : [nsd] x [nnd]
	 * \param det determinant of the transformation: [nip] */
	void ComputeDNa(const LocalArrayT& coords, ArrayT<dArray2DT>& DNa,
		dArrayT& det);

	/** compute nodal values.
	 * project the integration point values to the nodes 
	 * \param ipvalues field values from a single integration pt: [numvals]
	 * \param nodalvalues extrapolated values: [nnd] x [numvals]
	 * \param IPnum integration point number */
	void NodalValues(const dArrayT& IPvalues, dArray2DT& nodalvalues,
		int IPnum) const; 	

	/** print the shape function values to the output stream */
	void Print(ostream& out) const;
	
	/** number of element facets */
	int NumFacets(void) const;
	
	/** return the local node numbers per facet.
	 * local nodes are numbered to produce at outward normal in the 
	 * order: vertex nodes, mid-edge nodes, mid-face nodes
	 * \param facet element face number
	 * \param facetnodes local node number on the facet */
	void NodesOnFacet(int facet, iArrayT& facetnodes) const;

	/** return a list of the number of nodes on every element face */
	void NumNodesOnFacets(iArrayT& num_nodes) const;

	/** return the nodes on each facet needed to determine neighbors
	 * across facets */
	void NeighborNodeMap(iArray2DT& facetnodes) const;
	
	/** return geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geom, iArrayT& facet_nodes) const;
	
  private:

	/* dimensions */
	GeometryT::CodeT fGeometryCode; // geometry shape code
	int fNumSD;        // number of spatial dimensions
	int fNumIP;        // number of integration points
	int fNumNodes;     // number of domain nodes

	/* parent domain shape functions and derivatives */
	dArray2DT		  fNa;
	ArrayT<dArray2DT> fDNa;
	dArrayT           fWeights; // ip weights

	/* parent domain geometry */
	GeometryBaseT* fGeometry;

	/* work space */
	dMatrixT fNodalExtrap; // extrapolation matrix
	dMatrixT fJacobian;    // jacobian matrix
};

/* inlines */

/* number of spatial dimensions */
inline int ParentDomainT::NumSD(void)        const { return fNumSD;        }
inline int ParentDomainT::NumIP(void)        const { return fNumIP;        }
inline int ParentDomainT::NumNodes(void)     const { return fNumNodes;     }
inline GeometryT::CodeT ParentDomainT::GeometryCode(void) const { return fGeometryCode; }

/* access to domain shape functions */
inline const dArray2DT& ParentDomainT::Na(void) const { return fNa; }
inline const double* ParentDomainT::Shape(int IPnum) const { return fNa(IPnum); }
inline const double* ParentDomainT::DShape(int IPnum, int dim) const
{
	return (fDNa[IPnum])(dim);	
}
inline const double* ParentDomainT::Weight(void)     const { return fWeights.Pointer(); }

/* returns jacobian of the nodal values with respect
* to the variables of the shape function derivaties */
inline void ParentDomainT::DomainJacobian(const LocalArrayT& nodal, int numIP,
dMatrixT& jacobian) const
{
	Jacobian(nodal, fDNa[numIP] , jacobian);
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
inline int ParentDomainT::NumFacets(void) const
{
	return fGeometry->NumFacets();
}

inline void ParentDomainT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	fGeometry->NodesOnFacet(facet, facetnodes);
}

inline void ParentDomainT::NumNodesOnFacets(iArrayT& num_nodes) const
{
	fGeometry->NumNodesOnFacets(num_nodes);
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
inline void ParentDomainT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	fGeometry->NeighborNodeMap(facetnodes);
}

/* return geometry and number of nodes on each facet */
inline void ParentDomainT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geom,
	iArrayT& facet_nodes) const
{
	fGeometry->FacetGeometry(facet_geom, facet_nodes);
}

#endif /* _PARENT_DOMAIN_T_H_ */
