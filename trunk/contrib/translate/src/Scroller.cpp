/* $Id: Scroller.cpp,v 1.2 2004-11-16 01:04:35 paklein Exp $ */
#include "Scroller.h"
#include "ExceptionT.h"
#include "OutputSetT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

Scroller::Scroller(ostream& out, istream& in, bool write):
	TranslateIOManager(out, in, write),
	fOutputID(-1),
	fCleavagePlane(0.0),
	fDirection(0)
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
	
	/* dimensions */
	int nnv = fModel.NumNodeVariables();
	int nev = fModel.NumElementVariables();
	int nnd = fModel.NumNodes();
	int nel = fModel.NumElements();
	
	/* set changing coordinates */
	iArrayT map(nnd);
	fModel.AllNodeMap(map);
	dArray2DT coordinates = fModel.Coordinates();
	fOutput->SetCoordinates(coordinates, &map);
	LocalArrayT curr_coords(LocalArrayT::kUnspecified, 0, coordinates.MinorDim());
	curr_coords.SetGlobal(coordinates);

	/* original connectivities */
	const ArrayT<StringT>& ids = fModel.ElementGroupIDs();
	ArrayT<const iArray2DT*> connectivities;
	fModel.ElementGroupPointers(ids, connectivities);
	iArrayT element_nodes;

	/* work space */
	dArray2DT n_values(nnd, nnv), n_values_last;
	dArray2DT e_values, e_values_all(nel, nev);
	nVariArray2DT<double> e_values_man(0, e_values, nev);

//TEMP
if (fDirection != 1)
	ExceptionT::GeneralFail(caller, "propagation direction must be +1: %d", fDirection);

	/* run through increments */
	int incr = time_steps.Length()/10;
	incr = (incr < 1) ? 1 : incr;
	AutoArrayT<int> keep_block(nel,0), keep_all(nel,0);
	for (int j = 0; j < time_steps.Length(); j++) 
	{
		/* progress */
		if ((j+1)%incr == 0)
			cout << "step " << j+1 << endl;
	
	    /* read node values (across all blocks) */
		input->ReadAllNodeVariables(j, n_values);

	    /* read element values (across all blocks) */
		input->ReadAllElementVariables(j, e_values_all);

		/* check for scrolling */
		if (j > 0) {
			double x_shift = fXLeft - fMeshSize;
			bool shift = false;
			for (int i = 0; i < fNodes.Length(); i++) 
				if (fabs(n_values_last(fNodes[i], dy_index)) > fOpenLB && /* was "open" */
					fabs(n_values(fNodes[i], dy_index)) < fCloseUB && /* now "closed" */
				    coordinates(fNodes[i],0) > x_shift) /* point is farthest right */
				{
					x_shift = coordinates(fNodes[i],0);	
					shift = true;
				}
					
			
			/* shift coordinates */
			if (shift)
			{
				/* include one more column of elements */
				x_shift += fMeshSize;

				/* reset reference coordinates */
				double divider = fXLeft - fMeshSize/10.0;
				for (int i = 0; i < coordinates.MajorDim(); i++)
				{
					double& x = coordinates(i,0);
					x -= x_shift; /* scroll mesh */
					if (x < divider)
						x += fPeriodicLength; /* move ahead of crack */
				}
			}
		}
		n_values_last = n_values;
		
		/* remove the "back" element */
		int index_all = 0;
		keep_all.Dimension(0);
		for (int i = 0; i < connectivities.Length(); i++)
		{
			keep_block.Dimension(0);
			const iArray2DT& connects = *(connectivities[i]);
			int nen = connects.MinorDim();
			curr_coords.Dimension(nen, coordinates.MinorDim());
			for (int k = 0; k < connects.MajorDim(); k++)
			{
				/* collect current reference coordinates */
				connects.RowAlias(k, element_nodes);
				curr_coords.SetLocal(element_nodes);

				/* find element bounds */
				const double* px = curr_coords(0);
				double x_min = *px;
				double x_max = *px;
				px++;
				for (int l = 1; l < nen; l++) {
					/* bounds */
					x_min = (*px < x_min) ? *px : x_min;
					x_max = (*px > x_max) ? *px : x_max;
					px++;
				}
				
				/* not connecting ends of the mesh */
				if (x_max - x_min < fPeriodicLength/2.0) {
					keep_all.Append(index_all);
					keep_block.Append(k);
				}
				index_all++;
			}
			fConnectivities_man[i].SetMajorDimension(keep_block.Length(), false);
			fConnectivities[i].RowCollect(keep_block, connects);
		}
		e_values_man.SetMajorDimension(keep_all.Length(), false);
		e_values.RowCollect(keep_all, e_values_all);

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
	cout << "\n Source file: ";
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

	cout << "\n smallest displacement considered \"open\": ";
	fIn >> fOpenLB;
	if (fEcho) fEchoOut << fOpenLB << "\n";

	cout << "\n largest displacement considered \"closed\": ";
	fIn >> fCloseUB;
	if (fEcho) fEchoOut << fCloseUB << "\n";
	
	if (fOpenLB < fCloseUB)
		ExceptionT::GeneralFail(caller, "open lower bound < closed upper bound: %g < %g",
			fOpenLB, fCloseUB);

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

	/* get bounds and periodic distance */
	fMeshSize = 0.0;
	double x0 = coordinates(fNodes[0],0);
	fXLeft = fXRight = x0;
	for (int i = 1; i < fNodes.Length(); i++)
	{
		double x = coordinates(fNodes[i],0);
		double dx = fabs(x - x0);

		/* bounds */
		if (x < fXLeft)
			fXLeft = x;
		else if (x > fXRight)
			fXRight = x;
		
		/* mesh size */	
		if (dx > kSmall && (dx < fMeshSize || fMeshSize < kSmall))
			fMeshSize = dx;
	}
	fPeriodicLength = fXRight - fXLeft + fMeshSize;
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
	
	/* collect pointers to changing connectivities */
	int neg = connectivities.Length();
	fConnectivities.Dimension(neg);
	fConnectivities_man.Dimension(neg);
	ArrayT<const iArray2DT*> connectivities_tmp(neg);
	for (int i = 0; i < neg; i++) {
		const iArray2DT& connects = *(connectivities[i]);
		fConnectivities_man[i].SetWard(0, fConnectivities[i], connects.MinorDim());
		fConnectivities_man[i].SetMajorDimension(connects.MajorDim(), false);
		fConnectivities[i] = connects;
		connectivities_tmp[i] = fConnectivities.Pointer(i);
	}
	
	/* define output */
	bool changing_geometry = true;
	OutputSetT output_set(fModel.ElementGroupGeometry(ids[0]), ids, connectivities_tmp, 
		n_labels, e_labels, changing_geometry);

	/* register output */
	fOutputID = fOutput->AddElementSet(output_set);

	/* clean up */
	delete input;
}
