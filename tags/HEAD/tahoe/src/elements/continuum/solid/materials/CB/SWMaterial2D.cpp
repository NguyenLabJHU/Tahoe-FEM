/* $Id: SWMaterial2D.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (08/25/1996)                                          */

#include "SWMaterial2D.h"

/* constructor */
SWMaterial2D::SWMaterial2D(ifstreamT& in, const ElasticT& element):
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
