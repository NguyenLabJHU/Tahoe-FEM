/* $Id: DomainIntegrationT.h,v 1.3 2001-07-11 01:03:30 paklein Exp $ */
/* created: paklein (09/04/1998) */

#ifndef _DOMAIN_INTEGRATION_T_H_
#define _DOMAIN_INTEGRATION_T_H_

/* direct members */
#include "ParentDomainT.h" // needed for inlines and geometry codes
#include "iArrayT.h"

/** class to manage the parent domain. Includes construction for
 * shared parent domains, integration point iterations, and some
 * basic access to integration domain information. "copy" constructor
 * creates "linked" objects which (i) share the same parent domain and
 * (ii) are synchronized in integration through the current integration
 * point reference. */
class DomainIntegrationT
{
public:

	/** constructor. 
	 * \param geometry_code geometry of the parent domain
	 * \param numIP number of integration points 
	 * \param numnodes number of domain nodes */
	DomainIntegrationT(GeometryT::CodeT geometry_code, int numIP, int numnodes);

	/** constructor. 
	 * \param link shared parent domain and "synch-ed" CurrIP */
	DomainIntegrationT(const DomainIntegrationT& link);

	/** destructor */
	virtual ~DomainIntegrationT(void);

	/** class-dependent initializations */
	virtual void Initialize(void);

	/** weights for all the integration points */
	const double* IPWeights(void) const;

	/* accessors */
	int NumSD(void) const; /**< number of spatial dimensions */
	int NumIP(void) const; /**< number of integration points */

	/* integration control */
	void TopIP(void);   /**< restart loop over integration points */
	int  NextIP(void);  /**< next integration point. \return 0 when after last ip */
	void SetIP(int ip); /**< move to specified integration point */
	const int& CurrIP(void) const; /**< reference to the "current" integration point number */

	/** array nodal shape functions at the "current" integration point */
	const double* IPShape(void) const;

	/** integration weight of the "current" integration point */
	double IPWeight(void) const;

	/** extrapolate values from the "current" integration point to the nodes.
	 * \param IPvalues values from the integration point: [nval] 
	 * \param nodalvalues extrapolated values: [nnd] x [nval] */
	void Extrapolate(const dArrayT& IPvalues, dArray2DT& nodalvalues) const;

	/** print shape functions and derivatives to out */
	virtual void Print(ostream& out) const;

	/** return the number of domain facets */
	int  NumFacets(void) const;

	/** list of number of nodes on each domain facet */ 
	void NumNodesOnFacets(iArrayT& num_nodes) const;
	
	/** local node numbering over facets.
	 * \param facet cannonical domain facet number
	 * \param facetnodes list of number of nodes on each facet 
	 * \note facetnodes does not need to be dimensioned */
	void NodesOnFacet(int facet, iArrayT& facetnodes) const;

	/** geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geom,
		iArrayT& facet_nodes) const;

	/** the nodes on each facet needed to determine neighbors
	 * across facets. \note this list is generally shorter
	 * than the lists returned by NodesOnFacet */
	void NeighborNodeMap(iArray2DT& facetnodes) const;

	/** shape functions for the specified face */
	const ParentDomainT& FacetShapeFunction(int facet) const;
	
	/** reference to the parent domain */
	const ParentDomainT& ParentDomain(void) const;

protected:

	/** access to domain shape functions */
	const dArray2DT& Na(void) const;

private:

	/** set surface shapefunctions */
	void SetSurfaceShapes(void);

protected:

	/* integration management */
	int  fNumIP;  /**< number of integration points */
	int& fCurrIP; /**< current integration point number */
	
	/** parent geometry */
	ParentDomainT* fDomain;
	
	/** face shapefunctions */
	ArrayT<ParentDomainT*> fSurfShapes;

	/** flags to make duplicated face shape functions */
	iArrayT fDelete;

private:

	/* flags used to "link" domain to another */
	int fDeleteDomain;
	int frefCurrIP;
};

/* inlines */

/* data for all integration points at once */
inline const double* DomainIntegrationT::IPWeights(void) const
{
return fDomain->Weight();
}

/* accessors */
inline int DomainIntegrationT::NumSD(void) const { return fDomain->NumSD(); }
inline int DomainIntegrationT::NumIP(void) const { return fNumIP; }

/* integration management */
inline void DomainIntegrationT::TopIP(void)  { fCurrIP = -1; }
inline int  DomainIntegrationT::NextIP(void) { return (++fCurrIP < fNumIP); }
inline const int& DomainIntegrationT::CurrIP(void) const { return fCurrIP; }
inline void DomainIntegrationT::SetIP(int ip)
{
#if __option(extended_errorcheck)
	/* check */
	if (ip < 0 || ip >= fNumIP)
	{
		cout << "\n DomainIntegrationT::SetIP: " << ip
		     << " is out of range {" << 0 << ", " << fNumIP << "}" << endl;
		throw eOutOfRange;
	}
#endif
	fCurrIP = ip;
}

/* data for the current integration point */
inline const double* DomainIntegrationT::IPShape(void) const
{
	return fDomain->Shape(fCurrIP);
}

inline double DomainIntegrationT::IPWeight(void) const
{
#if __option(extended_errorcheck)
	/* range checking */
	if (fCurrIP < 0 || fCurrIP >= fNumIP) throw eOutOfRange;
#endif
	return *(fDomain->Weight() + fCurrIP);
}

/* extrapolate integration point values */
inline void DomainIntegrationT::Extrapolate(const dArrayT& IPvalues,
	dArray2DT& nodalvalues) const
{
	fDomain->NodalValues(IPvalues, nodalvalues, CurrIP());
}	

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
inline int DomainIntegrationT::NumFacets(void) const
{
	return fDomain->NumFacets();
}

inline void DomainIntegrationT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	fDomain->NodesOnFacet(facet, facetnodes);
}

inline void DomainIntegrationT::NumNodesOnFacets(iArrayT& num_nodes) const
{
	fDomain->NumNodesOnFacets(num_nodes);
}

/* return geometry and number of nodes on each facet */
inline void DomainIntegrationT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geom,
	iArrayT& facet_nodes) const
{
	fDomain->FacetGeometry(facet_geom, facet_nodes);
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
inline void DomainIntegrationT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	fDomain->NeighborNodeMap(facetnodes);
}

/* return shapefunctions for the specified element facet */
inline const ParentDomainT& DomainIntegrationT::FacetShapeFunction(int facet) const
{
	return *fSurfShapes[facet];
}

/* reference to the parent domain */
inline const ParentDomainT& DomainIntegrationT::ParentDomain(void) const
{
	if (!fDomain) throw eGeneralFail;
	return *fDomain;
}

/* access to domain shape functions */
inline const dArray2DT& DomainIntegrationT::Na(void) const
{	
	return fDomain->Na();
}

#endif /* _DOMAIN_INTEGRATION_T_H_ */