/* $Id: C1DomainIntegrationT.cpp,v 1.2 2003-10-08 21:04:55 rdorgan Exp $ */
/* class to manage the parent domain including construction for           */
/* shared parent domains, integration point iterations, and some          */
/* basic access to integration domain information. "copy" constructor     */
/* creates "linked" objects which (i) share the same parent domain and    */
/* (ii) are synchronized in integration through the current integration   */
/* point reference.                                                       */

#include "C1DomainIntegrationT.h"

/* constructor */

using namespace Tahoe;

C1DomainIntegrationT::C1DomainIntegrationT(C1GeometryT::CodeT geometry_code, int numIP, int numnodes):
        fNumIP(numIP),
        fCurrIP(frefCurrIP),
        fDeleteDomain(1),
        frefCurrIP(0)
{
        /* set parent geometry */
        fC1Domain = new C1ParentDomainT(geometry_code, fNumIP, numnodes);
        if (!fC1Domain) throw ExceptionT::kOutOfMemory;
        
        /* set parent domain shape functions and derivatives */
        fC1Domain->Initialize();

        /* set surface shapefunctions */
        if (fC1Domain->C1GeometryCode() != C1GeometryT::kC1Line) SetSurfaceShapes();
}

C1DomainIntegrationT::C1DomainIntegrationT(const C1DomainIntegrationT& link):
        fNumIP(link.fNumIP),
        fCurrIP(link.fCurrIP),
        fC1Domain(link.fC1Domain),
        fDeleteDomain(0),
        frefCurrIP(-1) // won't be used
{
        /* set surface shapefunctions */
        if (fC1Domain->C1GeometryCode() != C1GeometryT::kC1Line)
        {
                /* surface shapefunctions */
                fSurfShapes.Alias(link.fSurfShapes);

                fDelete.Dimension(fSurfShapes.Length());
                fDelete = 0;
        }
}

/* destructor */
C1DomainIntegrationT::~C1DomainIntegrationT(void)
{
        if (fDeleteDomain) delete fC1Domain;

        /* surface shapefunctions */
        for (int i = 0; i < fSurfShapes.Length(); i++)
                if (fDelete[i] == 1)
                {
                        delete fSurfShapes[i];
                        fSurfShapes[i] = NULL;
                }
}

/* set all local parameters */
void C1DomainIntegrationT::Initialize(void) {} //nothing to do

/* print the shape function values to the output stream */
void C1DomainIntegrationT::Print(ostream& out) const
{
        fC1Domain->Print(out);
}

/**********************************************************************
* Protected
**********************************************************************/

/* set surface shapefunctions */
void C1DomainIntegrationT::SetSurfaceShapes(void)
{
        /* memory */
        int num_facets = NumFacets();
        fSurfShapes.Dimension(num_facets);
        fDelete.Dimension(num_facets);
        
        /* surface shape information */
        ArrayT<C1GeometryT::CodeT> facet_geom;
        iArrayT num_facet_nodes;
        fC1Domain->FacetGeometry(facet_geom, num_facet_nodes);

        /* construct surface shapes */
        for (int i = 0; i < num_facets; i++)
        {
                C1GeometryT::CodeT geo = facet_geom[i];
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
                        fSurfShapes[i] = new C1ParentDomainT(geo, nip, nnd);
                        if (!fSurfShapes[i]) throw ExceptionT::kOutOfMemory;
        
                        /* set shape functions and derivatives */
                        fSurfShapes[i]->Initialize();

                        /* mark for deletion */
                        fDelete[i] = 1;
                }
        }
}
