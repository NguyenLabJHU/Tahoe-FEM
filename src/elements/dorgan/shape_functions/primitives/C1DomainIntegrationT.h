/* $Id: C1DomainIntegrationT.h,v 1.1 2003-09-29 19:59:03 rdorgan Exp $ */
#ifndef _C1_DOMAIN_INTEGRATION_T_H_
#define _C1_DOMAIN_INTEGRATION_T_H_

/* direct members */
#include "C1ParentDomainT.h" // needed for inlines and geometry codes
#include "iArrayT.h"

namespace Tahoe {

/** class to manage the parent domain. Includes construction for
 * shared parent domains, integration point iterations, and some
 * basic access to integration domain information. "copy" constructor
 * creates "linked" objects which (i) share the same parent domain and
 * (ii) are synchronized in integration through the current integration
 * point reference. */
class C1DomainIntegrationT
{
public:

	/** constructor. 
	 * \param geometry_code geometry of the parent domain
	 * \param numIP number of integration points 
	 * \param numnodes number of domain nodes */
	C1DomainIntegrationT(C1GeometryT::CodeT geometry_code, int numIP, int numnodes);

	/** constructor. 
	 * \param link shared parent domain and "synch-ed" CurrIP */
	C1DomainIntegrationT(const C1DomainIntegrationT& link);

	/** destructor */
	virtual ~C1DomainIntegrationT(void);

	/** class-dependent initializations */
	virtual void Initialize(void);

	/** weights for all the integration points */
	const double* IPWeights(void) const;

	/** \name accessors */
	/*@{*/
	int NumSD(void) const; /**< number of spatial dimensions */
	int NumIP(void) const; /**< number of integration points */
	C1GeometryT::CodeT C1GeometryCode(void) const; /**< domain geometry */
	/*@}*/

	/** \name integration control */
	/*@{*/
	void TopIP(void);   /**< restart loop over integration points */
	int  NextIP(void);  /**< next integration point. \return 0 when after last ip */
	void SetIP(int ip); /**< move to specified integration point */
	const int& CurrIP(void) const; /**< reference to the "current" integration point number */
	/*@}*/

	/** array nodal shape functions at the "current" integration point */
	const double* IPShape(void) const;

	/** integration weight of the "current" integration point */
	double IPWeight(void) const;

	/** extrapolate values from the "current" integration point to the nodes.
	 * \param IPvalues values from the integration point: [nval] 
	 * \param nodalvalues extrapolated values: [nnd] x [nval] */
	void Extrapolate(const dArrayT& IPvalues, dArray2DT& nodalvalues) const;

	/** extrapolate values the integration point values to the nodes.
	 * \param IPvalues values from the integration points: [nip] 
	 * \param nodalvalues extrapolated values: [nnd] */
	void ExtrapolateAll(const dArrayT& IPvalues, dArrayT& nodalvalues) const;

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
	void FacetGeometry(ArrayT<C1GeometryT::CodeT>& facet_geom,
		iArrayT& facet_nodes) const;

	/** the nodes on each facet needed to determine neighbors
	 * across facets. \note this list is generally shorter
	 * than the lists returned by NodesOnFacet */
	void NeighborNodeMap(iArray2DT& facetnodes) const;

	/** shape functions for the specified face */
	const C1ParentDomainT& FacetShapeFunction(int facet) const;
	
	/** reference to the parent domain */
	const C1ParentDomainT& C1ParentDomain(void) const;

	/** evaluate the shape functions and gradients. Compute the values of the
	 * shape functions and their gradients at an arbirary point in the
	 * in the parent domain. Coordinates must fall within the domain.
	 * \param coords point in the parent domain
	 * \param Na destination for shape function values for each of the domain
	 *        nodes. Must be dimensioned: [nnd]
	 * \param DNa destination for shape function derivatives. Must be 
	 *        dimensioned: [nsd] x [nnd] */
	void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, 
		dArray2DT& DNa) const;

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
	C1ParentDomainT* fC1Domain;
	
	/** face shapefunctions */
	ArrayT<C1ParentDomainT*> fSurfShapes;

	/** flags to make duplicated face shape functions */
	iArrayT fDelete;

private:

	/* flags used to "link" domain to another */
	int fDeleteDomain;
	int frefCurrIP;
};

/* inlines */

/* data for all integration points at once */
inline const double* C1DomainIntegrationT::IPWeights(void) const
{
return fC1Domain->Weight();
}

/* accessors */
inline int C1DomainIntegrationT::NumSD(void) const { return fC1Domain->NumSD(); }
inline int C1DomainIntegrationT::NumIP(void) const { return fNumIP; }
inline C1GeometryT::CodeT C1DomainIntegrationT::C1GeometryCode(void) const 
{ 
	return fC1Domain->C1GeometryCode(); 
}

/* integration management */
inline void C1DomainIntegrationT::TopIP(void)  { fCurrIP = -1; }
inline int  C1DomainIntegrationT::NextIP(void) { return (++fCurrIP < fNumIP); }
inline const int& C1DomainIntegrationT::CurrIP(void) const { return fCurrIP; }
inline void C1DomainIntegrationT::SetIP(int ip)
{
#if __option(extended_errorcheck)
	/* check */
	if (ip < 0 || ip >= fNumIP)
	{
		cout << "\n C1DomainIntegrationT::SetIP: " << ip
		     << " is out of range {" << 0 << ", " << fNumIP << "}" << endl;
		throw ExceptionT::kOutOfRange;
	}
#endif
	fCurrIP = ip;
}

/* data for the current integration point */
inline const double* C1DomainIntegrationT::IPShape(void) const
{
	return fC1Domain->Shape(fCurrIP);
}

inline double C1DomainIntegrationT::IPWeight(void) const
{
#if __option(extended_errorcheck)
	/* range checking */
	if (fCurrIP < 0 || fCurrIP >= fNumIP) throw ExceptionT::kOutOfRange;
#endif
	return *(fC1Domain->Weight() + fCurrIP);
}

/* extrapolate integration point values */
inline void C1DomainIntegrationT::Extrapolate(const dArrayT& IPvalues,
	dArray2DT& nodalvalues) const
{
	fC1Domain->NodalValues(IPvalues, nodalvalues, CurrIP());
}	

inline void C1DomainIntegrationT::ExtrapolateAll(const dArrayT& IPvalues,
	dArrayT& nodalvalues) const
{
	fC1Domain->NodalValues(IPvalues, nodalvalues);
}	

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
inline int C1DomainIntegrationT::NumFacets(void) const
{
	return fC1Domain->NumFacets();
}

inline void C1DomainIntegrationT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	fC1Domain->NodesOnFacet(facet, facetnodes);
}

inline void C1DomainIntegrationT::NumNodesOnFacets(iArrayT& num_nodes) const
{
	fC1Domain->NumNodesOnFacets(num_nodes);
}

/* return geometry and number of nodes on each facet */
inline void C1DomainIntegrationT::FacetGeometry(ArrayT<C1GeometryT::CodeT>& facet_geom,
	iArrayT& facet_nodes) const
{
	fC1Domain->FacetGeometry(facet_geom, facet_nodes);
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
inline void C1DomainIntegrationT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	fC1Domain->NeighborNodeMap(facetnodes);
}

/* return shapefunctions for the specified element facet */
inline const C1ParentDomainT& C1DomainIntegrationT::FacetShapeFunction(int facet) const
{
	return *fSurfShapes[facet];
}

/* reference to the parent domain */
inline const C1ParentDomainT& C1DomainIntegrationT::C1ParentDomain(void) const
{
	if (!fC1Domain) throw ExceptionT::kGeneralFail;
	return *fC1Domain;
}

/* evaluate the shape functions and gradients. */
inline void C1DomainIntegrationT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, 
	dArray2DT& DNa) const
{
	fC1Domain->EvaluateShapeFunctions(coords, Na, DNa);
}

/* access to domain shape functions */
inline const dArray2DT& C1DomainIntegrationT::Na(void) const
{	
	return fC1Domain->Na();
}

} // namespace Tahoe 
#endif /* _C1_DOMAIN_INTEGRATION_T_H_ */
