/* $Id: SWMaterial2D.cpp,v 1.5.30.1 2004-01-21 19:10:02 paklein Exp $ */
/* created: paklein (08/25/1996) */
#include "SWMaterial2D.h"

using namespace Tahoe;

/* constructor */
SWMaterial2D::SWMaterial2D(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("Stillinger_Weber_material_2D"),
	NL_E_RotMat2DT(in, support, kPlaneStrain),
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
