/* $Id: SWMaterial2D.cpp,v 1.2 2001-07-03 01:35:00 paklein Exp $ */
/* created: paklein (08/25/1996)                                          */

#include "SWMaterial2D.h"

/* constructor */
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
