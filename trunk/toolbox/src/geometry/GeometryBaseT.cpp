/* $Id: GeometryBaseT.cpp,v 1.4 2003-11-10 22:14:29 cjkimme Exp $ */
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
