/* $Id: C1GeometryBaseT.cpp,v 1.3 2003-11-19 20:38:24 rdorgan Exp $ */
#include "C1GeometryBaseT.h"
#include "ExceptionT.h"
#include <iostream.h>

using namespace Tahoe;

/* constructor */
C1GeometryBaseT::C1GeometryBaseT(int numnodes): //, int numfacets):
	fNumNodes(numnodes)
        //	fNumFacets(numfacets)
{


}

/* destructor */
C1GeometryBaseT::~C1GeometryBaseT(void) { }

        ///* compute gradients of the "bubble" modes */
        //void C1GeometryBaseT::BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const
        //{
        //#pragma unused(Na_x)
        //        cout << "\n C1GeometryBaseT::BubbleModeGradients: geometry does not have bubble modes\n" 
        //             <<   "     Derived classes must override to define bubble mode derivatives" 
        //             << endl;
        //        throw ExceptionT::kGeneralFail;
        //}
