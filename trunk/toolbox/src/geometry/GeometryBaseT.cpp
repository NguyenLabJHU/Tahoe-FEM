/* $Id: GeometryBaseT.cpp,v 1.8 2005-01-26 19:52:10 paklein Exp $ */
/* created: paklein (10/21/1997) */
#include "GeometryBaseT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
GeometryBaseT::GeometryBaseT(int numnodes, int numfacets):
	fNumNodes(numnodes),
	fNumFacets(numfacets)
{

}

/* destructor */
GeometryBaseT::~GeometryBaseT(void) { }

/* compute gradients of the "bubble" modes */
void GeometryBaseT::BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const
{
#pragma unused(Na_x)
	ExceptionT::GeneralFail("GeometryBaseT::BubbleModeGradients", 
		"no bubble modes for geometry \"%s\"", ToString(Geometry()));
}

/* return true if the given point is within the domain defined by */
bool GeometryBaseT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
#pragma unused(coords)
#pragma unused(point)
	ExceptionT::GeneralFail("GeometryBaseT::PointInDomain", 
		"not implemented for geometry \"%s\"", ToString(Geometry()));
	return false;
}

/* integration point gradient matrix */
void GeometryBaseT::IPGradientTransform(int ip, dMatrixT& transform) const
{
#pragma unused (ip)
#pragma unused (transform)
	ExceptionT::GeneralFail("GeometryBaseT::IPGradientTransform", 
		"not implemented for geometry \"%s\"", ToString(Geometry()));
}

/* return the integration point whose domain contains the given point in the
 * parent domain coordinates */
int GeometryBaseT::IPDomain(int nip, const dArrayT& coords) const
{
#pragma unused (nip)
#pragma unused (coords)
	ExceptionT::GeneralFail("GeometryBaseT::IPDomain", 
		"not implemented for geometry \"%s\"", ToString(Geometry()));
	return -1;
}

/* subdomain geometry */
GeometryT::CodeT GeometryBaseT::NodalSubDomainGeometry(void) const {
	ExceptionT::GeneralFail("GeometryBaseT::NodalSubDomainGeometry", 
		"not implemented for geometry \"%s\"", ToString(Geometry()));
	return GeometryT::kNone;
}

/* number of nodes defining the nodal subdomain */
int GeometryBaseT::NodalSubDomainNumPoints(void) const {
	ExceptionT::GeneralFail("GeometryBaseT::NodalSubDomainNumPoints", 
		"not implemented for geometry \"%s\"", ToString(Geometry()));
	return -1;
}
	
/* compute the coordinates of the points defining the nodal subdomain */
void GeometryBaseT::NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
	LocalArrayT& subdomain_coords) const
{
#pragma unused(coords)
#pragma unused(node)
#pragma unused(subdomain_coords)
	ExceptionT::GeneralFail("GeometryBaseT::NodalSubDomainCoordinates", 
		"not implemented for geometry \"%s\"", ToString(Geometry()));
}
