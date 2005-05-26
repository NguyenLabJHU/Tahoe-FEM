/* $Id: DomainIntegrationT.cpp,v 1.5 2004-02-02 23:45:00 paklein Exp $ */
/* created: paklein (09/04/1998)                                          */
/* class to manage the parent domain including construction for           */
/* shared parent domains, integration point iterations, and some          */
/* basic access to integration domain information. "copy" constructor     */
/* creates "linked" objects which (i) share the same parent domain and    */
/* (ii) are synchronized in integration through the current integration   */
/* point reference.                                                       */

#include "DomainIntegrationT.h"

/* constructor */

using namespace Tahoe;

DomainIntegrationT::DomainIntegrationT(GeometryT::CodeT geometry_code, int numIP, int numnodes):
	fNumIP(numIP),
	fCurrIP(frefCurrIP),
	fDeleteDomain(1),
	frefCurrIP(-1)
{
	/* set parent geometry */
	fDomain = new ParentDomainT(geometry_code, fNumIP, numnodes);
	if (!fDomain) throw ExceptionT::kOutOfMemory;
	
	/* set parent domain shape functions and derivatives */
	fDomain->Initialize();

	/* set surface shapefunctions */
	if (fDomain->GeometryCode() != GeometryT::kLine) SetSurfaceShapes();
}

DomainIntegrationT::DomainIntegrationT(const DomainIntegrationT& link):
	fNumIP(link.fNumIP),
	fCurrIP(link.fCurrIP),
	fDomain(link.fDomain),
	fDeleteDomain(0),
	frefCurrIP(-1) // won't be used
{
	/* set surface shapefunctions */
	if (fDomain->GeometryCode() != GeometryT::kLine)
	{
		/* surface shapefunctions */
		fSurfShapes.Alias(link.fSurfShapes);

		fDelete.Dimension(fSurfShapes.Length());
		fDelete = 0;
	}
}

/* destructor */
DomainIntegrationT::~DomainIntegrationT(void)
{
	if (fDeleteDomain) delete fDomain;

	/* surface shapefunctions */
	for (int i = 0; i < fSurfShapes.Length(); i++)
		if (fDelete[i] == 1)
		{
			delete fSurfShapes[i];
			fSurfShapes[i] = NULL;
		}
}

/* set all local parameters */
void DomainIntegrationT::Initialize(void) {} //nothing to do

/* print the shape function values to the output stream */
void DomainIntegrationT::Print(ostream& out) const
{
	fDomain->Print(out);
}

/**********************************************************************
* Protected
**********************************************************************/

/* set surface shapefunctions */
void DomainIntegrationT::SetSurfaceShapes(void)
{
	/* memory */
	int num_facets = NumFacets();
	fSurfShapes.Dimension(num_facets);
	fDelete.Dimension(num_facets);
	
	/* surface shape information */
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT num_facet_nodes;
	fDomain->FacetGeometry(facet_geom, num_facet_nodes);

	/* construct surface shapes */
	for (int i = 0; i < num_facets; i++)
	{
		GeometryT::CodeT geo = facet_geom[i];
		int nnd = num_facet_nodes[i];
	
		int same_shape = -1;
		for (int j = 0; j < i && same_shape < 0; j++)
			if (facet_geom[j] == geo &&
			    num_facet_nodes[j] == nnd) same_shape = j;	

		/* duplicate surface shape */
		if (same_shape > -1)
		{
			fSurfShapes[i] = fSurfShapes[same_shape];
			fDelete[i] = 0;
		}
		/* construct new surface shape */
		else
		{
			//TEMP - set number of ip's on facet same as number of nodes
			int nip = (NumSD() == 3) ? 4 : nnd;
			fSurfShapes[i] = new ParentDomainT(geo, nip, nnd);
			if (!fSurfShapes[i]) throw ExceptionT::kOutOfMemory;
	
			/* set shape functions and derivatives */
			fSurfShapes[i]->Initialize();

			/* mark for deletion */
			fDelete[i] = 1;
		}
	}
}
