/* $Id: SWMaterial2D.cpp,v 1.2.6.1 2002-06-27 18:03:01 cjkimme Exp $ */
/* created: paklein (08/25/1996)                                          */

#include "SWMaterial2D.h"

/* constructor */

using namespace Tahoe;

SWMaterial2D::SWMaterial2D(ifstreamT& in, const FiniteStrainT& element):
	NL_E_RotMat2DT(in, element, kPlaneStrain),
	SWDataT(in)
{

}

/* I/O functions */
void SWMaterial2D::Print(ostream& out) const
{
	/* inherited */
	NL_E_RotMat2DT::Print(out);
	SWDataT::Write(out); //SWDataT
}
