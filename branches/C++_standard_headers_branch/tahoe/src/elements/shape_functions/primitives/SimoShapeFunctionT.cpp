/* $Id: SimoShapeFunctionT.cpp,v 1.8 2005-07-14 07:11:07 paklein Exp $ */

#include "SimoShapeFunctionT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* constructor */
SimoShapeFunctionT::SimoShapeFunctionT(GeometryT::CodeT geometry_code, 
	int numIP, const LocalArrayT& coords, const LocalArrayT& element_modes):
	ShapeFunctionT(geometry_code, numIP, coords),
	fElementModes(element_modes),
	fHas3DIncompressibleMode(false)
{
	/* check coordinates type */
	if (fCoords.Type() != LocalArrayT::kInitCoords)
	{
		cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting local reference coordinates: "
		     << fCoords.Type() << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* checks */
	if (geometry_code == GeometryT::kQuadrilateral)
	{
		/* number of element nodes */
		if (coords.NumberOfNodes() != 4)
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 4 element nodes: "
			     << coords.NumberOfNodes() << endl;
			throw ExceptionT::kBadInputValue;
		}
		
		/* check number of enhanced modes */
		if (element_modes.NumberOfNodes() != 2)
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 2 enhanced modes: "
			     << element_modes.NumberOfNodes() << endl;
			throw ExceptionT::kBadInputValue;
		}
	}
	else if (geometry_code == GeometryT::kHexahedron)
	{
		/* number of element nodes */
		if (coords.NumberOfNodes() != 8)
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 4 element nodes: "
			     << coords.NumberOfNodes() << endl;
			throw ExceptionT::kBadInputValue;
		}

		/* check number of enhanced modes */
		if (element_modes.NumberOfNodes() != 3) //TEMP - no incompressible mode yet
		{
			cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: expecting 3 enhanced modes: "
			     << element_modes.NumberOfNodes() << endl;
			throw ExceptionT::kBadInputValue;
			
			/* set flag */
			fHas3DIncompressibleMode = true; //TEMP - only set if num_modes = 4
		}
	}
	else
	{
		cout << "\n SimoShapeFunctionT::SimoShapeFunctionT: element geometry must be\n" 
		     <<   "     quad/hex for 2D/3D: " << geometry_code << endl;
		throw ExceptionT::kBadInputValue;
	}

	/* allocate derivatives of bubble modes */
	int num_bubble_modes = NumSD();
	int num_derivatives  = NumSD();
	fDNaX_bubble.Dimension(NumIP());
	fDNa_bubble.Dimension(NumIP());
	for (int i = 0; i < NumIP(); i++)
	{
		fDNaX_bubble[i].Dimension(num_derivatives, num_bubble_modes);
		fDNa_bubble[i].Dimension(num_derivatives, num_bubble_modes);
	}
	
	/* 3D incompressible mode */	
	if (fHas3DIncompressibleMode)
	{
		/* derivatives of 1 mode */
		fDNaX_inc.Dimension(NumIP(), num_derivatives);
	}
	
	/* allocate work space */
	fNa_0.Dimension(fCoords.NumberOfNodes());
	fDNa_0.Dimension(NumSD(), fCoords.NumberOfNodes());
	fJ.Dimension(NumSD());
	fJ_0_inv.Dimension(NumSD());
}

/* initialization class. */
void SimoShapeFunctionT::Initialize(void)
{
	/* inherited */
	ShapeFunctionT::Initialize();
	
	/* reference to the parent domain geometry */
	const GeometryBaseT& geometry = ParentDomain().Geometry();
	
	/* compute/store gradients of the bubble modes */
	geometry.BubbleModeGradients(fDNa_bubble);
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
		throw ExceptionT::kGeneralFail;
	}
	else
	{
		/* inherited - set Galerkin part */
		ShapeFunctionT::SetDerivatives();

		/* number of integration points */
		int nip = NumIP();
		
		/* for 2D:5/3D:9 point rules, assume the final point is at the origin */
		bool has_center_point = (NumSD() == 2 && nip == 5) || (NumSD() == 3 && nip == 9);
	
		/* Jacobian at the origin of the parent domain */
		if (has_center_point)
		{
			/* assume last ip is at the origin */
			fDomain->DomainJacobian(fCoords, nip - 1, fJ_0_inv);
		}
		else
		{
			/* compute shape function derivatives at the origin */
			double px[3] = {0.0, 0.0, 0.0};
			dArrayT coords_0(NumSD(), px);
			fDomain->EvaluateShapeFunctions(coords_0, fNa_0, fDNa_0);

			/* compute jacobian */
			fDomain->Jacobian(fCoords, fDNa_0, fJ_0_inv);
		}
		
		/* set Jacobians */
		double j_0 = fJ_0_inv.Det();
		fJ_0_inv.Inverse();
		
		/* transform derivatives */
		int i_centroid = (has_center_point) ? nip - 1 : -1;
		for (int i = 0; i < nip; i++)
		{
			/* apply change of variables to the shape function derivatives */
			TransformDerivatives(fJ_0_inv, fDNa_bubble[i], fDNaX_bubble[i]);
			
			/* skip for point at centroid */
			if (i != i_centroid)
			{
				/* compute jacobian */
				fDomain->DomainJacobian(fCoords, i, fJ);

				/* scale */
				fDNaX_bubble[i] *= j_0/fJ.Det();
			}
		}
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
