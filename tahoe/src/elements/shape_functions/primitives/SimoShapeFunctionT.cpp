/* $Id: SimoShapeFunctionT.cpp,v 1.1 2001-07-11 01:03:30 paklein Exp $ */
/* created: paklein (06/26/1996)                                          */

#include "SimoShapeFunctionT.h"
#include "LocalArrayT.h"

/* constructor */
SimoShapeFunctionT::SimoShapeFunctionT(GeometryT::CodeT geometry_code, 
	int numIP, const LocalArrayT& coords, const dArray2DT& element_modes):
	ShapeFunctionT(geometry_code, numIP, coords, kStandardB),
	fElementModes(element_modes)
{
//number of nodes should be 4 for 2D/8 for 3D
//number of integration points should be 4/5 for 2D or
//8/9 for 3D
//number of modes should be 2 x 2 in 2D or
//3/4 x 3 in 3D
}

/* compute local shape functions and derivatives */ 	
void SimoShapeFunctionT::SetDerivatives(void)
{

// compute all shape function gradients

}

/* print the shape function values to the output stream */
void SimoShapeFunctionT::Print(ostream& out) const
{
	/* inherited */
	ShapeFunctionT::Print(out);

	out << "\n Enhanced bubble mode shape function derivatives:\n";
	for (int i = 0; i < fDNaX_bubble.Length(); i++)
		fDNaX_bubble[i].WriteNumbered(out);
		
	if (fDNaX_inc.Length() > 0)
	{
		out << "\n Enhanced incompressible mode shape function derivatives:\n";
		for (int i = 0; i < fDNaX_inc.Length(); i++)
			fDNaX_inc[i].WriteNumbered(out);
	}	
}
