/* $Id: CSEBaseT.cpp,v 1.6.6.1 2001-10-26 15:27:51 sawimme Exp $ */
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
const int CSEBaseT::NumNodalOutputCodes = 5;
const int CSEBaseT::NumElementOutputCodes = 3;

/* constructor */
CSEBaseT::CSEBaseT(FEManagerT& fe_manager):
	ElementBaseT(fe_manager),
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
	fNodalOutputCodes.Allocate(NumNodalOutputCodes);
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();
	for (int i = 0; i < fNodalOutputCodes.Length(); i++)
	{
		in >> fNodalOutputCodes[i];

		/* convert all to "at print increment" */
		if (fNodalOutputCodes[i] != IOBaseT::kAtNever)
			fNodalOutputCodes[i] = IOBaseT::kAtInc;		
	}

	/* echo */
	out << " Number of nodal output codes. . . . . . . . . . = " << fNodalOutputCodes.Length() << '\n';
	out << "    [" << fNodalOutputCodes[NodalCoord   ] << "]: initial nodal coordinates\n";
	out << "    [" << fNodalOutputCodes[NodalDisp    ] << "]: nodal displacements\n";
	out << "    [" << fNodalOutputCodes[NodalDispJump] << "]: nodal gap (magnitude)\n";
	out << "    [" << fNodalOutputCodes[NodalTraction] << "]: nodal traction magnitude\n";
	out << "    [" << fNodalOutputCodes[MaterialData ] << "]: constitutive output data\n";

	/* check */
	if (fNodalOutputCodes.Min() < IOBaseT::kAtFinal ||
	    fNodalOutputCodes.Max() > IOBaseT::kAtInc) throw eBadInputValue;

	fElementOutputCodes.Allocate(NumElementOutputCodes);
	fElementOutputCodes = IOBaseT::kAtNever;

//TEMP - backward compatibility
	if (StringT::versioncmp(fFEManager.Version(), "v3.01") < 1)
	{
		/* message */
		cout << "\n CSEBaseT::Initialize: use input file version newer than v3.01\n" 
		     <<   "     to enable element output control" << endl;
		out << "\n CSEBaseT::Initialize: use input file version newer than v3.01\n" 
		    <<   "     to enable element output control" << endl;
	}
	else
	{
//TEMP - BACK
		int num_codes = fElementOutputCodes.Length();
		if (StringT::versioncmp(fFEManager.Version(), "v3.02") < 1)
		{
			cout << "\n CSEBaseT::Initialize: use input file version newer than v3.02\n"
		         <<   "     to enable output control of element averaged traction" << endl;
			out << "\n CSEBaseT::Initialize: use input file version newer than v3.02\n"
		         <<   "     to enable output control of element averaged traction" << endl;
			num_codes = 2;
		}
	
		/* read in at a time to allow comments */
		for (int j = 0; j < num_codes; j++)
		{
			in >> fElementOutputCodes[j];
		
			/* convert all to "at print increment" */
			if (fElementOutputCodes[j] != IOBaseT::kAtNever)
				fElementOutputCodes[j] = IOBaseT::kAtInc;
		}	

		/* checks */
		if (fElementOutputCodes.Min() < IOBaseT::kAtFail ||
		    fElementOutputCodes.Max() > IOBaseT::kAtInc) throw eBadInputValue;
	}

	/* echo */
	out << " Number of element output codes. . . . . . . . . = " << fElementOutputCodes.Length() << '\n';
	out << "    [" << fElementOutputCodes[Centroid      ] << "]: centroid\n";
	out << "    [" << fElementOutputCodes[CohesiveEnergy] << "]: dissipated cohesive energy\n";
	out << "    [" << fElementOutputCodes[Traction      ] << "]: average traction\n";
	
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

	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);
	iArrayT block_ID(fBlockData.MajorDim());
	fBlockData.ColumnCopy(kID, block_ID);

	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

	/* set output specifier */
	int ID = fFEManager.ElementGroupNumber(this) + 1;
	OutputSetT output_set(ID, geo_code, block_ID, fConnectivities, n_labels, 
		e_labels, false);
		
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

	/* map output flags to count of values */
	iArrayT n_counts;
	SetNodalOutputCodes(mode, fNodalOutputCodes, n_counts);
	iArrayT e_counts;
	SetElementOutputCodes(mode, fElementOutputCodes, e_counts);

	/* calculate output values */
	dArray2DT n_values;
	dArray2DT e_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);

	/* send to output */
	fFEManager.WriteOutput(fOutputID, n_values, e_values);
}

/* compute specified output parameter and send for smoothing */
void CSEBaseT::SendOutput(int kincode)
{
	/* output flags */
	iArrayT flags(fNodalOutputCodes.Length());

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
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

	/* reset averaging workspace */
	fNodes->ResetAverage(n_counts.Sum());

	/* set flags for no element output */
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	/* generate nodal values */
	dArray2DT e_values, n_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
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

void CSEBaseT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
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

void CSEBaseT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Allocate(flags.Length());
	counts = 0;

	if (flags[Centroid] == mode)
		counts[Centroid] = fNumSD;
	if (flags[CohesiveEnergy] == mode)
		counts[CohesiveEnergy] = 1;
	if (flags[Traction] == mode)
		counts[Traction] = 1;
}

/* construct output labels array */
void CSEBaseT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	/* allocate nodal output labels */
	n_labels.Allocate(n_codes.Sum());

	int count = 0;
	if (n_codes[NodalDisp])
	{
		if (fNumDOF > 3) throw eGeneralFail;
		const char* dlabels[3] = {"D_X", "D_Y", "D_Z"};

		for (int i = 0; i < fNumDOF; i++)
			n_labels[count++] = dlabels[i];
	}

	if (n_codes[NodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < fNumSD; i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[NodalDispJump]) n_labels[count++] = "jump";
	if (n_codes[NodalTraction]) n_labels[count++] = "Tmag";
	
	/* allocate nodal output labels */
	e_labels.Allocate(e_codes.Sum());
	count = 0;
	if (e_codes[Centroid])
	{
		const char* xlabels[] = {"xc_1", "xc_2", "xc_3"};
		for (int i = 0; i < fNumSD; i++)
			e_labels[count++] = xlabels[i];
	}
	if (e_codes[CohesiveEnergy]) e_labels[count++] = "phi";
	if (e_codes[Traction]) e_labels[count++] = "Tmag";
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
	        const iArrayT& elemnodes = fElementCards[i].NodesX();
		int* nodes = elemnodes.Pointer();

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
