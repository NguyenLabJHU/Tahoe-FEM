/* $Id: Scroller.cpp,v 1.1 2004-11-15 09:21:19 paklein Exp $ */
#include "Scroller.h"
#include "ExceptionT.h"
#include "OutputSetT.h"

using namespace Tahoe;

Scroller::Scroller(ostream& out, istream& in, bool write):
	TranslateIOManager(out, in, write),
	fOutputID(-1),
	fCleavagePlane(0.0),
	fDirection(0),
	fJumpThreshold(0.0)
{

}

void Scroller::Translate(const StringT& program, const StringT& version, const StringT& title)
{
	const char caller[] = "Scroller::Translate";

	/* set sources */
	SetInput();
	SetOutput(program, version, title);

	/* write geometry */
	WriteGeometry();

	/* define output sets */	
	InitOutputSets();

	/* source */
	InputBaseT* input = IOBaseT::NewInput(fFileType, cout);
	input->Open(fSource);
	
	/* get time increments */
	dArrayT time_steps(input->NumTimeSteps());
	input->ReadTimeSteps(time_steps);

	/* look for y-displacement in output */
	ArrayT<StringT> nlabels;
	input->ReadNodeLabels(nlabels);
	int dy_index = -1;
	for (int j = 0; j < nlabels.Length(); j++)
		if (nlabels[j] == "D_Y")
			dy_index = j;
	if (dy_index == -1)
		ExceptionT::GeneralFail(caller, "\"D_Y\" not found in nodal output");

	/* set changing coordinates */
	iArrayT map;
	fModel.AllNodeMap(map);
	dArray2DT coordinates = fModel.Coordinates();
	fOutput->SetCoordinates(coordinates, &map);

	/* work space */
	int nnv = fModel.NumNodeVariables();
	int nev = fModel.NumElementVariables();
	int nnd = fModel.NumNodes();
	int nel = fModel.NumElements();
	dArray2DT n_values(nnd, nnv), n_values_last;
	dArray2DT e_values(nel, nev);

//TEMP
ExceptionT::GeneralFail(caller, "propagation direction must be +1: %d", fDirection);

	/* run through increments */
	for (int j = 0; j < time_steps.Length(); j++) 
	{
	    /* read node values (across all blocks) */
		input->ReadAllNodeVariables(j, n_values);

	    /* read element values (across all blocks) */
		input->ReadAllElementVariables(j, e_values);

		/* check for scrolling */
		if (j > 0) {
			int hit = -1;
			for (int i = 0; i < fNodes.Length(); i++)
			{
				double jump = n_values(fNodes[i], dy_index) - n_values_last(fNodes[i], dy_index);
				if (fabs(jump) > fJumpThreshold)
				{
					//check for hit == -1
					//keep track of right most jump
				}
			}
			
			//process hit
		}
		n_values_last = n_values;

		/* write to output */
		fOutput->WriteOutput(time_steps[j], fOutputID, n_values, e_values);
	}

	/* clean up */
	delete input;
}

/************************************************************************
 * Protected
 ************************************************************************/

void Scroller::SetInput(void)
{
	const char caller[] = "Scroller::SetInput";

	/* source file */
	cout << "\n Source files: ";
	fIn >> fSource;
	if (fEcho) fEchoOut << fSource << "\n";

	/* tracking information */
	cout << "\n y-coordinate of cleavage plane: ";
	fIn >> fCleavagePlane;
	if (fEcho) fEchoOut << fCleavagePlane << "\n";

	int count = 0;
	while (++count < 5 && fDirection != 1 && fDirection != -1) {
		cout << "\n propagation direction (-1: left, +1: right): ";
		fIn >> fDirection;
	}
	if (fDirection != 1 && fDirection != -1)
		ExceptionT::GeneralFail(caller, "bad direction %d", fDirection);
	else if (fEcho) 
		fEchoOut << fDirection << "\n";
	
	cout << "\n displacement jump threshold: ";
	fIn >> fJumpThreshold;
	if (fEcho) fEchoOut << fJumpThreshold << "\n";
	
	/* try to guess file format */
	fFileType = IOBaseT::name_to_FileTypeT(fSource);

	/* open the first file to transfer the model geometry */
	fSource.ToNativePathName();
	if (!fModel.Initialize(fFileType, fSource, true))
		ExceptionT::DatabaseFail(caller, "unable to initialize model file %s", 
			fSource.Pointer());

	/* look for nodes on the cleavage plane */
	AutoArrayT<int> nodes(fModel.NumNodes()/10, 25);
	nodes.Dimension(0);
	const dArray2DT& coordinates = fModel.Coordinates();
	for (int i = 0; i < coordinates.MajorDim(); i++)
		if (fabs(coordinates(i,1) - fCleavagePlane) < kSmall)
			nodes.Append(i);

	if (nodes.Length() == 0)
		ExceptionT::GeneralFail(caller, "no nodes found on the cleavage plane y = %g", fCleavagePlane);
	fNodes.Dimension(nodes.Length());
	nodes.CopyInto(fNodes);
}

/************************************************************************
 * Private
 ************************************************************************/

/* init output sets */
void Scroller::InitOutputSets(void)
{
	const char caller[] = "Scroller::InitOutputSets";
	if (!fOutput) ExceptionT::GeneralFail(caller);

	/* use first source to define output set */
	InputBaseT* input = IOBaseT::NewInput(fFileType, cout);
	input->Open(fSource);
	
	/* construct new output set */
	const ArrayT<StringT>& ids = fModel.ElementGroupIDs();
	ArrayT<const iArray2DT*> connectivities;
	fModel.ElementGroupPointers(ids, connectivities);
	ArrayT<StringT> n_labels, e_labels;
	fModel.NodeLabels(n_labels);
	fModel.ElementLabels(e_labels);
	
	/* define output */
	bool changing_geometry = true;
	OutputSetT output_set(fModel.ElementGroupGeometry(ids[0]), ids, connectivities, 
		n_labels, e_labels, changing_geometry);

	/* register output */
	fOutputID = fOutput->AddElementSet(output_set);

	/* clean up */
	delete input;
}
