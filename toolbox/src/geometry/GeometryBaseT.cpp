/* $Id: GeometryBaseT.cpp,v 1.7 2004-05-17 05:11:04 paklein Exp $ */
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
