/* $Id: CSEBaseT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (11/19/1997)                                          */

#include "CSEBaseT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "SurfaceShapeT.h"
#include "NodeManagerT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"

/* initialize static data */
const int CSEBaseT::NumOutputCodes = 5;

/* constructor */
CSEBaseT::CSEBaseT(FEManagerT& fe_manager):
	ElementBaseT(fe_manager),
	fOutputCodes(NumOutputCodes),
	fLocInitCoords1(LocalArrayT::kInitCoords),
	fLocCurrCoords(LocalArrayT::kCurrCoords),
	fFractureArea(0.0),
	fShapes(NULL)
{
	/* read control parameters */
	ifstreamT& in = fFEManager.Input();

	in >> fGeometryCode;
	in >> fNumIntPts;
	in >> fCloseSurfaces;
	in >> fOutputArea;

	/* checks */
	if (fNumSD == 2 && fGeometryCode != GeometryT::kLine)
	{
		cout << "\n CSEBaseT::CSEBaseT: expecting geometry code "
		     << GeometryT::kLine<< " for 2D: " << fGeometryCode << endl;
		throw eBadInputValue;
	}
	else if (fNumSD == 3 &&
	         fGeometryCode != GeometryT::kQuadrilateral &&
	         fGeometryCode != GeometryT::kTriangle)
	{
		cout << "\n CSEBaseT::CSEBaseT: expecting geometry code " << GeometryT::kQuadrilateral
		     << " or\n" <<   "     " << GeometryT::kTriangle << " for 3D: "
		     << fGeometryCode << endl;
		throw eBadInputValue;
	}
	
	if (fCloseSurfaces != 0 &&
	    fCloseSurfaces != 1) throw eBadInputValue;
	if (fOutputArea != 0 &&
	    fOutputArea != 1) throw eBadInputValue;
}

/* destructor */
CSEBaseT::~CSEBaseT(void)
{
	delete fShapes;
	fShapes = NULL;
}

/* allocates space and reads connectivity data */
void CSEBaseT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* dimensions */
	int num_facet_nodes = fNumElemNodes/2;

	/* initialize local arrays */
	fLocInitCoords1.Allocate(num_facet_nodes, fNumSD);
	fLocCurrCoords.Allocate(fNumElemNodes, fNumSD);
	fFEManager.RegisterLocal(fLocInitCoords1);
	fFEManager.RegisterLocal(fLocCurrCoords);

	/* construct surface shape functions */
	fShapes = new SurfaceShapeT(fGeometryCode, fNumIntPts, fNumElemNodes, fNumDOF,
		fLocInitCoords1);
	if (!fShapes) throw eOutOfMemory;
	fShapes->Initialize();

	/* work space */
	fNodes1.Allocate(num_facet_nodes);
	fNEEvec.Allocate(fNumElemEqnos);
	fNEEmat.Allocate(fNumElemEqnos);

	/* echo output codes (one at a time to allow comments) */
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();
	for (int i = 0; i < fOutputCodes.Length(); i++)
	{
		in >> fOutputCodes[i];

		/* convert all to "at print increment" */
		if (fOutputCodes[i] != IOBaseT::kAtNever)
			fOutputCodes[i] = IOBaseT::kAtInc;		
	}

	out << " Number of nodal output codes. . . . . . . . . . = " << NumOutputCodes << '\n';
	out << "    [" << fOutputCodes[NodalCoord   ] << "]: initial nodal coordinates\n";
	out << "    [" << fOutputCodes[NodalDisp    ] << "]: nodal displacements\n";
	out << "    [" << fOutputCodes[NodalDispJump] << "]: nodal gap (magnitude)\n";
	out << "    [" << fOutputCodes[NodalTraction] << "]: nodal traction magnitude\n";
	out << "    [" << fOutputCodes[MaterialData ] << "]: constitutive output data\n";

	/* check */
	if (fOutputCodes.Min() < IOBaseT::kAtFinal ||
	    fOutputCodes.Max() > IOBaseT::kAtInc) throw eBadInputValue;
	
	/* close surfaces */
	if (fCloseSurfaces) CloseSurfaces();

	/* output stream */
	if (fOutputArea == 1)
	{
		/* generate file name */
		StringT name = (fFEManager.Input()).filename();
		name.Root();
		name.Append(".grp", fFEManager.ElementGroupNumber(this) + 1);
		name.Append(".fracture");
		
		/* open stream */
		farea_out.open(name);
	}
}

/* initial condition/restart functions (per time sequence) */
void CSEBaseT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();

	/* initialize element status flags */
	for (int i = 0; i < fNumElements; i++)
		fElementCards[i].Flag() = kON;
}

/* finalize time increment */
void CSEBaseT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* deactivate marked elements */
	for (int i = 0; i < fNumElements; i++)
	{
		int& flag = fElementCards[i].Flag();
		flag = (flag == kMarked) ? kOFF : flag;
	}
}

/* resets to the last converged solution */
void CSEBaseT::ResetStep(void)
{
	/* inherited */
	ElementBaseT::ResetStep();

	/* unset marks */
	for (int i = 0; i < fNumElements; i++)
	{
		int& flag = fElementCards[i].Flag();
		flag = (flag == kMarked) ? kON : flag;
	}
}

/* solution calls */
void CSEBaseT::AddNodalForce(int node, dArrayT& force)
{
//TEMP
#pragma unused(node)
#pragma unused(force)
	//not implemented
}

/* returns the energy as defined by the derived class types */
double CSEBaseT::InternalEnergy(void) { return 0.0; } //not implemented

/* writing output */
void CSEBaseT::RegisterOutput(void)
{
//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"

	/* ID is just group number */
	int ID = fFEManager.ElementGroupNumber(this) + 1;

	/* "deformed" geometry */
	GeometryT::CodeT geo_code;
	switch (fGeometryCode)
{
		case GeometryT::kLine:		
			geo_code = GeometryT::kQuadrilateral;
			break;

	case GeometryT::kQuadrilateral:
			geo_code = GeometryT::kHexahedron;
			break;

		case GeometryT::kTriangle:
			geo_code = GeometryT::kPentahedron;
			break;

		default:	

		cout << "\n CSEBaseT::RegisterOutput: could not translate\n";
		cout << "     geometry code " << fGeometryCode
			 << " to a pseudo-geometry code for the volume." << endl;
		throw eGeneralFail;	
	}	

	/* get output configuration */
	iArrayT counts;
	SetOutputCodes(IOBaseT::kAtInc, fOutputCodes, counts);
	int num_out = counts.Sum();

	/* variable labels */
	ArrayT<StringT> n_labels(num_out);
	ArrayT<StringT> e_labels;
	GenerateOutputLabels(counts, n_labels);

	OutputSetT output_set(ID, geo_code, fConnectivities,
		n_labels, e_labels, false);
		
	/* register and get output ID */
	fOutputID = fFEManager.RegisterOutput(output_set);
}

//NOTE - this function is identical to ContinuumElementT::WriteOutput
void CSEBaseT::WriteOutput(IOBaseT::OutputModeT mode)
{
//TEMP - not handling general output modes yet
	if (mode != IOBaseT::kAtInc)
	{
		cout << "\n CSEBaseT::WriteOutput: only handling \"at increment\"\n"
		     <<   "     print mode. SKIPPING." << endl;
		return;
	}

	/* fracture area */
	if (fOutputArea)
	{
		farea_out << setw(kDoubleWidth) << fFEManager.Time();
		farea_out << setw(kDoubleWidth) << fFractureArea << endl;
	}

	/* nodal output */
	iArrayT codes;
	SetOutputCodes(mode, fOutputCodes, codes);
	int num_out = codes.Sum();

	dArray2DT group_n_values;
	dArray2DT group_e_values(0,0);
	if (num_out > 0)
	{
		/* reset averaging workspace */
		fNodes->ResetAverage(num_out);

		/* compute nodal values */
		ComputeNodalValues(codes);

		/* get nodal values */
		fNodes->OutputUsedAverage(group_n_values);		
	}

	/* send out */
	fFEManager.WriteOutput(fOutputID, group_n_values, group_e_values);

}

/* compute specified output parameter and send for smoothing */
void CSEBaseT::SendOutput(int kincode)
{
/* output flags */
iArrayT flags(fOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case NodalDisp:
		    flags[NodalDisp] = fNumDOF;
			break;
		case NodalDispJump:
		    flags[NodalDispJump] = 1;
			break;
		case NodalTraction:
		    flags[NodalTraction] = 1;
			break;
		default:
			cout << "\n CSEBaseT::SendKinematic: invalid output code: ";
			cout << kincode << endl;
	}

	/* number of output values */
	iArrayT counts;
	SetOutputCodes(IOBaseT::kAtInc, flags, counts);
	int num_out = counts.Sum();

	/* reset averaging workspace */
	fNodes->ResetAverage(num_out);

	/* generate nodal values */
	ComputeNodalValues(counts);
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void CSEBaseT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	/* control parameters */
	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode << '\n';
	out << "    eq." << GeometryT::kLine          << ", line\n";
	out << "    eq." << GeometryT::kQuadrilateral << ", quadrilateral\n";
	out << "    eq." << GeometryT::kTriangle	  << ", triangle\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIntPts     << '\n';
	out << " Initial surface closure flag. . . . . . . . . . = " << fCloseSurfaces << '\n';
	out << " Output fracture surface area. . . . . . . . . . = " << fOutputArea    << '\n';
}

void CSEBaseT::SetOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const
{
	/* initialize */
	counts.Allocate(flags.Length());
	counts = 0;

	if (flags[NodalCoord] == mode)
		counts[NodalCoord] = fNumSD;
	if (flags[NodalDisp] == mode)
		counts[NodalDisp] = fNumDOF;
	if (flags[NodalDispJump] == mode)
		counts[NodalDispJump] = 1;
	if (flags[NodalTraction] == mode)
		counts[NodalTraction] = 1;
}

/* write all current element information to the stream */
void CSEBaseT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;

	out <<   " initial coords:\n";
	temp.Allocate(fLocInitCoords1.NumberOfNodes(), fLocInitCoords1.MinorDim());
	fLocInitCoords1.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " current coords:\n";
	temp.Allocate(fLocCurrCoords.NumberOfNodes(), fLocCurrCoords.MinorDim());
	fLocCurrCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* construct output labels array */
void CSEBaseT::GenerateOutputLabels(const iArrayT& codes,
	ArrayT<StringT>& labels) const
{
	/* allocate */
	labels.Allocate(codes.Sum());

	int count = 0;
	if (codes[NodalDisp])
	{
		if (fNumDOF > 3) throw eGeneralFail;
		const char* dlabels[3] = {"D_X", "D_Y", "D_Z"};

		for (int i = 0; i < fNumDOF; i++)
			labels[count++] = dlabels[i];
	}

	if (codes[NodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};

		for (int i = 0; i < fNumSD; i++)
			labels[count++] = xlabels[i];
	}

	if (codes[NodalDispJump]) labels[count++] = "jump";
	if (codes[NodalTraction]) labels[count++] = "Tmag";
}

/* close surfaces to zero gap */
void CSEBaseT::CloseSurfaces(void) const
{
	/* get coordinates */
	const dArray2DT& init_coords = fNodes->InitialCoordinates();

	/* local nodes numbers on each facet */
	const iArray2DT& facetnodes = fShapes->NodesOnFacets();
		
	/* collapse elements */
	for (int i = 0; i < fNumElements; i++)
	{			
		int* pfacet1 = facetnodes(0);
		int* pfacet2 = facetnodes(1);
		int*   nodes = fConnectivities(i);

		for (int j = 0; j < facetnodes.MinorDim(); j++)
		{
			/* facet coordinates */		
			double* px1 = init_coords(nodes[*pfacet1++]);
			double* px2 = init_coords(nodes[*pfacet2++]);
				
			for (int k = 0; k < fNumSD; k++)
			{
				double x_mid = 0.5*(*px1 + *px2);
				*px1++ = x_mid;
				*px2++ = x_mid;
			}
		}		
	}
}

/* return the default number of element nodes */
int CSEBaseT::DefaultNumElemNodes(void) const
{
	/* return number of element nodes given the facet geometry */
	switch (fGeometryCode)
	{
		case GeometryT::kLine:
			return 4;
		case GeometryT::kQuadrilateral:
			return 8;
		case GeometryT::kTriangle:
			return 6;
		default:
			cout << "\n CSEBaseT::DefaultNumElemNodes: unknown geometry code: "
			     << fGeometryCode << endl;
			return 0;
	}
}
//NOTE: needed because ExodusII does not store ANY information about
//      empty element groups, which causes trouble for parallel execution
//      when a partition contains no element from a group.
