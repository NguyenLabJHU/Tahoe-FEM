/* $Id: ParentDomainT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
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

class ParentDomainT
{
public:

	/* constructor */
	ParentDomainT(GeometryT::CodeT geometry_code, int numIP, int numnodes);

	/* destructor */
	~ParentDomainT(void);

	/* set all local parameters */
	void Initialize(void);

	/* accessors */
	int NumSD(void) const;
	int NumIP(void) const;
	int NumNodes(void) const;
	GeometryT::CodeT GeometryCode(void) const;

	/* access to domain shape functions */
	const dArray2DT& Na(void) const;      // all the shape function data
	const double* Shape(int IPnum) const; // nodal shape functions at the IP
	const double* DShape(int IPnum, int dim) const; // derivatives at the IP for the given dim
const double* Weight(void) const;     // integration weights for all IP's

	/* interpolation to the current integration point */
	void Interpolate(const LocalArrayT& nodal, dArrayT& interp, int IPnum) const;
	
	/* interpolate all to integration points: (nip x nu) */
	void Interpolate(const LocalArrayT& nodal, dArray2DT& interp) const;

	/* returns jacobian of the nodal values with respect
	 * to the variables of the shape function derivatives.
	 * Q returns as the transformation from global to local(')
	 * coordinates, i.e., t'_i = Q_ik t_k, where t'_j (j = nsd)
	 * is the "normal" direction */
	void Jacobian(const LocalArrayT& nodal, const dArray2DT& DNa,
		dMatrixT& jacobian) const;
	void DomainJacobian(const LocalArrayT& nodal, int numIP, dMatrixT& jacobian) const;
	double SurfaceJacobian(const dMatrixT& jacobian) const;
	double SurfaceJacobian(const dMatrixT& jacobian, dMatrixT& Q) const;

	/* chain rule jacobian of shape functions wrt coordinates that
	 * are passed in, for all integration points at once */
	void ComputeDNa(const LocalArrayT& coords, ArrayT<dArray2DT>& DNa,
		dArrayT& det);

	/* compute nodal values:
	 * ipvalues[numvals] : field values from a single integration pt
	 * nodalvalues[fNumNodes x numvals] : extrapolated values */
	void NodalValues(const dArrayT& IPvalues, dArray2DT& nodalvalues,
		int IPnum) const; 	

	/* print the shape function values to the output stream */
	void Print(ostream& out) const;
	
	/* return the local node numbers for each facet of the element
	 * numbered to produce at outward normal in the order: vertex
	 * nodes, mid-edge nodes, mid-face nodes */
	int  NumFacets(void) const;
	void NodesOnFacet(int facet, iArrayT& facetnodes) const;
	void NumNodesOnFacets(iArrayT& num_nodes) const;

	/* returns the nodes on each facet needed to determine neighbors
	 * across facets */
	void NeighborNodeMap(iArray2DT& facetnodes) const;
	
	/* return geometry and number of nodes on each facet */
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
