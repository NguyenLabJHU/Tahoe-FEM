/* $Id: DomainIntegrationT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (09/04/1998)                                          */
/* class to manage the parent domain including construction for           */
/* shared parent domains, integration point iterations, and some          */
/* basic access to integration domain information. "copy" constructor     */
/* creates "linked" objects which (i) share the same parent domain and    */
/* (ii) are synchronized in integration through the current integration   */
/* point reference.                                                       */

#ifndef _DOMAIN_INTEGRATION_T_H_
#define _DOMAIN_INTEGRATION_T_H_

/* direct members */
#include "ParentDomainT.h" // needed for inlines and geometry codes
#include "iArrayT.h"

class DomainIntegrationT
{
public:

	/* constructors */
	DomainIntegrationT(GeometryT::CodeT geometry_code, int numIP, int numnodes);
	DomainIntegrationT(const DomainIntegrationT& link);
		// synch-ed integration domains with shared ParentDomainT

	/* destructor */
	virtual ~DomainIntegrationT(void);

	/* class-dependent initializations */
	virtual void Initialize(void);

	/* weights for all the integration points */
const double* IPWeights(void) const;

	/* accessors */
	int NumSD(void) const;	
	int NumIP(void) const;

	/* integration management */
	void TopIP(void);
	int  NextIP(void);
	const int& CurrIP(void) const;
	void SetIP(int ip);

/**** data for the current integration point ****/
	const double* IPShape(void) const; // nodalshape functions
	double IPWeight(void) const;       // integration point weight
/************************************************/

	/* print shape functions and derivatives */
	virtual void Print(ostream& out) const;

	/* return the local node numbers for each facet of the element
	 * numbered to produce at outward normal in the order: vertex
	 * nodes, mid-edge nodes, mid-face nodes */
	int  NumFacets(void) const;
	void NodesOnFacet(int facet, iArrayT& facetnodes) const;
	void NumNodesOnFacets(iArrayT& num_nodes) const;

	/* return geometry and number of nodes on each facet */
	void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geom,
		iArrayT& facet_nodes) const;

	/* returns the nodes on each facet needed to determine neighbors
	 * across facets */
	void NeighborNodeMap(iArray2DT& facetnodes) const;

	/* return shapefunctions for the specified element facet */
	const ParentDomainT& FacetShapeFunction(int facet) const;
	
	/* reference to the parent domain */
	const ParentDomainT& ParentDomain(void) const;

protected:

	/* access to domain shape functions */
	const dArray2DT& Na(void) const;

private:

	/* set surface shapefunctions */
	void SetSurfaceShapes(void);

protected:

	/* integration management */
int  fNumIP;  // number of integration points
	int& fCurrIP; // current integration point
	
	/* parent geometry */
	ParentDomainT* fDomain;
	
	/* surface shapefunctions */
	ArrayT<ParentDomainT*> fSurfShapes;
	iArrayT fDelete;

private:

	int fDeleteDomain;
	int frefCurrIP; // for synching domains
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

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
inline int DomainIntegrationT::NumFacets(void) const
{
	return fDomain->NumFacets();
}

inline void DomainIntegrationT::NodesOnFacet(int facet,
	iArrayT& facetnodes) const
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
