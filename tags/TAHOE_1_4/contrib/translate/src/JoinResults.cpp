/* $Id: JoinResults.cpp,v 1.1 2004-05-10 01:28:48 paklein Exp $ */
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
			bool write_step = false;
			if (i == 0 ||(i > 0 && time_steps[j] > last_time)) {
				write_step = true;
				num_written++;
			}

			/* run through groups */
			for (int k = 0; write_step && k < ids.Length(); k++) {
			
			    /* read node values */
				input->ReadNodeVariables(j, ids[k], n_values);

			    /* read element values */
			    input->ReadElementVariables(j, ids[k], e_values);
		
				/* concat to output */
				fOutput->WriteOutput(time_steps[j], fOutputID[k], n_values, e_values);
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

	/* use first source to define output sets */
	InputBaseT* input = IOBaseT::NewInput(fFileTypes[0], cout);
	input->Open(fSources[0]);

	/* generate one output set per block ID */
	const ArrayT<StringT>& ids = fModel.ElementGroupIDs();
	fOutputID.Dimension(ids.Length());
	for (int i = 0; i < ids.Length(); i++)
	{
		/* construct new output set */
		ArrayT<StringT> block_ID(1);
		block_ID[0] = ids[i];
		ArrayT<const iArray2DT*> connectivities(1);
		connectivities[0] = fModel.ElementGroupPointer(ids[i]);
  		ArrayT<StringT> n_labels(fModel.NumNodeVariables());
		ArrayT<StringT> e_labels(fModel.NumElementVariables());
		fModel.NodeLabels(n_labels);
		fModel.ElementLabels(e_labels);
		
		/* define output */
		OutputSetT output_set(fModel.ElementGroupGeometry(ids[i]), block_ID, connectivities, 
			n_labels, e_labels, false);

		/* register output */
		fOutputID[i] = fOutput->AddElementSet(output_set);
	}
	
	/* clean up */
	delete input;
}
