/* $Id: GeometryBaseT.cpp,v 1.5 2004-02-28 21:52:26 paklein Exp $ */
/* created: paklein (10/21/1997) */
#include "GeometryBaseT.h"
#include "ExceptionT.h"
#include <iostream.h>

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
	ExceptionT::GeneralFail("GeometryBaseT::BubbleModeGradients", "no bubble modes for geometry %d", int(Geometry()));
}

/* return true if the given point is within the domain defined by */
bool GeometryBaseT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
#pragma unused(coords)
#pragma unused(point)
	ExceptionT::GeneralFail("GeometryBaseT::PointInDomain", "not implemented for geometry %d", int(Geometry()));
	return false;
}