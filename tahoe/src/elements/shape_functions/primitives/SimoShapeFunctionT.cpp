/* $Id: SimoShapeFunctionT.cpp,v 1.2 2001-07-13 02:19:10 paklein Exp $ */

#include "SimoShapeFunctionT.h"
#include "LocalArrayT.h"

/* constructor */
SimoShapeFunctionT::SimoShapeFunctionT(GeometryT::CodeT geometry_code, 
	int numIP, const LocalArrayT& coords, const dArray2DT& element_modes):
	ShapeFunctionT(geometry_code, numIP, coords, kStandardB),
	fElementModes(element_modes),
	fHas3DIncompressibleMode(false)
{
	/* checks */
	if (geometry_code == GeometryT::kQuadrilateral)
	{
		/* number of element nodes */
		if (coords.NumberOfNodes() != 4)
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 4 element nodes: "
			     << coords.NumberOfNodes() << endl;
			throw eBadInputValue;
		}
		
		/* check number of enhanced modes */
		if (element_modes.MajorDim() != 2)
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 2 enhanced modes: "
			     << element_modes.MajorDim() << endl;
			throw eBadInputValue;
		}
	}
	else if (geometry_code == GeometryT::kHexahedron)
	{
		/* number of element nodes */
		if (coords.NumberOfNodes() != 8)
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 4 element nodes: "
			     << coords.NumberOfNodes() << endl;
			throw eBadInputValue;
		}

		/* check number of enhanced modes */
		if (element_modes.MajorDim() != 3) //TEMP - no incompressible mode yet
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 3 enhanced modes: "
			     << element_modes.MajorDim() << endl;
			throw eBadInputValue;
			
			/* set flag */
			fHas3DIncompressibleMode = true; //TEMP - only set if num_modes = 4
		}
	}
	else
	{
		cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: element geometry must be\n" 
		     <<   "     quad/hex for 2D/3D: " << geometry_code << endl;
		throw eBadInputValue;
	}

	/* allocate derivatives of bubble modes */
	int num_bubble_modes = NumSD();
	int num_derivatives  = NumSD();
	fDNaX_bubble.Allocate(NumIP());
	for (int i = 0; i < NumIP(); i++)
		fDNaX_bubble[i].Allocate(num_derivatives, num_bubble_modes);
	
	/* 3D incompressible mode */	
	if (fHas3DIncompressibleMode)
	{
		/* derivatives of 1 mode */
		fDNaX_inc.Allocate(NumIP(), num_derivatives);
	}
}

/* compute local shape functions and derivatives */ 	
void SimoShapeFunctionT::SetDerivatives(void)
{
	/* 3D incompressible mode */	
	if (fHas3DIncompressibleMode)
	{
		//TEMP
		cout << "\n SimoShapeFunctionT::SetDerivatives: 3D incompressible mode\n"
		     <<   "    not supported yet" << endl;
		throw eGeneralFail;
	}
	else
	{
		/* inherited - set Galerkin part */
		ShapeFunctionT::SetDerivatives();
	
		//compute derivatives of enhanced modes in current config - here or will this
		//   happen within the shape functions themselves
		//compute enhancement to the deformation gradient
	}
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
		fDNaX_inc.WriteNumbered(out);
	}	
}
