/* $Id: GeometryBaseT.cpp,v 1.3 2002-10-20 22:32:08 paklein Exp $ */
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
	cout << "\n GeometryBaseT::BubbleModeGradients: geometry does not have bubble modes\n" 
	     <<   "     Derived classes must override to define bubble mode derivatives" 
	     << endl;
	throw ExceptionT::kGeneralFail;
}
