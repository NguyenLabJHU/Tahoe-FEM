/* $Id: GeometryBaseT.cpp,v 1.2 2002-09-30 20:52:43 sawimme Exp $ */
/* created: paklein (10/21/1997) */

#include "GeometryBaseT.h"
#include "ExceptionCodes.h"
#include <iostream.h>

/* constructor */

using namespace Tahoe;

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
	throw eGeneralFail;
}
