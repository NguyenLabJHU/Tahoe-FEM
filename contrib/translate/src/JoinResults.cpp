/* $Id: JoinResults.cpp,v 1.2 2004-11-11 03:57:24 paklein Exp $ */
#include "JoinResults.h"
#include "ExceptionT.h"
#include "OutputSetT.h"

using namespace Tahoe;

JoinResults::JoinResults(ostream& out, istream& in, bool write):
	TranslateIOManager(out, in, write)
{

}

void JoinResults::Translate(const StringT& program, const StringT& version, const StringT& title)
{
	const char caller[] = "JoinResults::Translate";

	/* set sources */
	SetInput();
	SetOutput(program, version, title);

	/* write geometry */
	WriteGeometry();

	/* define output sets */	
	InitOutputSets();

	/* collect time increments */
	double last_time = -99;
	const ArrayT<StringT>& ids = fModel.ElementGroupIDs();	
	for (int i = 0; i < fSources.Length(); i++) {
	
		/* source */
		InputBaseT* input = IOBaseT::NewInput(fFileTypes[i], cout);
		input->Open(fSources[i]);
	
		/* get time increments */
		dArrayT time_steps(input->NumTimeSteps());
		input->ReadTimeSteps(time_steps);

		/* dimensions */
		int nnv = fModel.NumNodeVariables();
		int nev = fModel.NumElementVariables();
		int nnd = fModel.NumNodes();
		int nel = fModel.NumElements();

		/* work space */
		dArray2DT n_values(nnd, nnv);
		dArray2DT e_values(nel, nev);

		/* run through increments */
		int num_written = 0;
		for (int j = 0; j < time_steps.Length(); j++) {

			/* steps must be sequential */
			if (i == 0 || (i > 0 && time_steps[j] > last_time)) {

			    /* read node values (across all blocks) */
				input->ReadAllNodeVariables(j, n_values);

			    /* read element values (across all blocks) */
				input->ReadAllElementVariables(j, e_values);

				/* concat to output */
				fOutput->WriteOutput(time_steps[j], fOutputID, n_values, e_values);
				num_written++;				
			}
		}
		
		/* report */
		cout << fSources[i] << ": " << num_written << " steps"<< endl;
		
		/* keep last time */
		if (time_steps.Length() > 0) last_time = time_steps.Last();

		/* clean up */
		delete input;
	}
}

/**************** PROTECTED **********************/

void JoinResults::SetInput(void)
{
	const char caller[] = "JoinResults::SetInput";

	int num_files = 0;
	cout << "\n Number of source files (> 1): ";
	fIn >> num_files;
	if (fEcho) fEchoOut << num_files << "\n";
	if (num_files < 2) ExceptionT::GeneralFail(caller, "Must have more than 1 file for merging.");
	
	fSources.Dimension(num_files);
	fFileTypes.Dimension(num_files);
	for (int i = 0; i < fSources.Length(); i++)
	{
		StringT database;
		cout << "file " << i+1 << ": ";
		fIn >> database;
		if (fEcho) fEchoOut << database << "\n";
		database.ToNativePathName();
	
		/* try to guess file format */
		fFileTypes[i] = IOBaseT::name_to_FileTypeT(database);
	
		/* store */
		fSources[i] = database;
    }    

	/* open the first file to transfer the model geometry */
	if (!fModel.Initialize(fFileTypes[0], fSources[0], true))
		ExceptionT::DatabaseFail(caller, "unable to initialize model file %s", 
			fSources[0].Pointer());
}

/**************** PRIVATE **********************/

/* init output sets */
void JoinResults::InitOutputSets(void)
{
	const char caller[] = "JoinResults::InitOutputSets";
	if (fSources.Length() < 1) ExceptionT::GeneralFail(caller);
	if (!fOutput) ExceptionT::GeneralFail(caller);

	/* use first source to define output set */
	InputBaseT* input = IOBaseT::NewInput(fFileTypes[0], cout);
	input->Open(fSources[0]);
	
	/* construct new output set */
	const ArrayT<StringT>& ids = fModel.ElementGroupIDs();
	ArrayT<const iArray2DT*> connectivities;
	fModel.ElementGroupPointers(ids, connectivities);
	ArrayT<StringT> n_labels, e_labels;
	fModel.NodeLabels(n_labels);
	fModel.ElementLabels(e_labels);
	
	/* define output */
	OutputSetT output_set(fModel.ElementGroupGeometry(ids[0]), ids, connectivities, 
		n_labels, e_labels, false);

	/* register output */
	fOutputID = fOutput->AddElementSet(output_set);

	/* clean up */
	delete input;
}
