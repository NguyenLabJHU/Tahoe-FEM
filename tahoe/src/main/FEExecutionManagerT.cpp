/* $Id: FEExecutionManagerT.cpp,v 1.69.2.9 2004-09-27 23:35:36 paklein Exp $ */
/* created: paklein (09/21/1997) */
#include "FEExecutionManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>
#include <ctype.h>
#include <stdlib.h>

#if defined(__MWERKS__) && __option(profile)
#include <Profiler.h>
#endif

/* element configuration header */
#include "ElementsConfig.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "FEManagerT.h"
#include "IOManager_mpi.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"
#include "StringT.h"
#include "GraphT.h"
#include "PartitionT.h"
#include "OutputSetT.h"
#include "ModelFileT.h"
#include "ExodusT.h"
#include "JoinOutputT.h"
#include "dArrayT.h"
#include "OutputBaseT.h"
#include "CommunicatorT.h"
#include "DecomposeT.h"

/* parameters */
#include "ParameterListT.h"
#include "ParameterTreeT.h"
#include "XML_Attribute_FormatterT.h"
#include "DotLine_FormatterT.h"

/* needed for bridging calculations FEExecutionManagerT::RunBridging */
#ifdef BRIDGING_ELEMENT
#include "FEManagerT_bridging.h"
#include "MultiManagerT.h"
#ifdef __DEVELOPMENT__
#include "BridgingScaleManagerT.h"
#include "FEManagerT_THK.h"
#endif
#include "TimeManagerT.h"
#include "NodeManagerT.h"
#include "dSPMatrixT.h"
#include "FieldT.h"
#include "IntegratorT.h"
#include "ElementBaseT.h"
#include "EAMFCC3D.h"
#endif /* BRIDGING_ELEMENT */

using namespace Tahoe;

/* Constructor */
FEExecutionManagerT::FEExecutionManagerT(int argc, char* argv[], char job_char,
	char batch_char, CommunicatorT& comm):
	ExecutionManagerT(argc, argv, job_char, batch_char, comm, 0)
{
#ifdef __CPLANT__
	/* if not prescribed as joined, write separate files */
	if (!CommandLineOption("-join_io")) AddCommandLineOption("-split_io");
#endif

	/* set communicator log level */
	if (CommandLineOption("-verbose")) comm.SetLogLevel(CommunicatorT::kLow);

	/* check for -dtd */
	if (CommandLineOption("-dtd")) AddCommandLineOption("-dtd");

	/* check for -xsd */
	if (CommandLineOption("-xsd")) AddCommandLineOption("-xsd");
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* MUST be overloaded */
void FEExecutionManagerT::RunJob(ifstreamT& in, ostream& status)
{
	const char caller[] = "FEExecutionManagerT::RunJob";
	
	/* set the path to the root file */
	StringT root;
	root.FilePath(in.filename());
	fstreamT::SetRoot(root);

	/* mode - job by default */
	ModeT mode = kJob;

	/* check command line options */
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
	{
		if (fCommandLineOptions[i] == "-decomp")
			mode = kDecompose;
		else if (fCommandLineOptions[i] == "-join")
			mode = kJoin;
	}

	/* check second char */
	if (mode == kJob) {
		
		/* peek at next char */
		char a = in.next_char();
		if (a == fJobChar) {
			mode = kBridging;
			in >> a; /* clear character */	
		}
	}

	switch (mode)
	{
		case kJob:
		{
			if (fComm.Size() == 1) {
				cout << "\n RunJob_analysis, serial version: " << in.filename() << endl;
				RunJob_analysis(in.filename(), status);
			} else {
				cout << "\n RunJob_analysis, parallel version: " << in.filename() << endl;
				RunJob_analysis(in.filename(), status);
			}
			break;
		}
		case kDecompose:
		{
			cout << "\n RunDecomp_serial: " << in.filename() << endl;
			if (fComm.Size() > 1) cout << " RunDecomp_serial: SERIAL ONLY" << endl;

			/* 'serial' communicator */
			int rank = fComm.Rank();
			CommunicatorT comm(fComm, (rank == 0) ? rank : CommunicatorT::kNoColor);
			
			/* decompose on rank 0 */
			if (rank == 0) RunDecomp_serial(in.filename(), status, comm);

			/* synch */
			fComm.Barrier();
			break;
		}
		case kJoin:
		{
			cout << "\n RunJoin_serial: " << in.filename() << endl;
			if (fComm.Size() > 1) cout << " RunJoin_serial: SERIAL ONLY" << endl;
			
			/* 'serial' communicator */
			int rank = fComm.Rank();
			CommunicatorT comm(fComm, (rank == 0) ? rank : CommunicatorT::kNoColor);

			/* join using rank 0 */
			if (rank == 0) RunJoin_serial(in.filename(), status);

			/* synch */
			fComm.Barrier();
			break;
		}
		default:
			ExceptionT::GeneralFail("FEExecutionManagerT::RunJob", "unknown mode: %d", mode);
	}

	/* clear the path to the root file */
	fstreamT::SetRoot(NULL);
}

/**********************************************************************
 * Protected
 **********************************************************************/

bool FEExecutionManagerT::AddCommandLineOption(const char* str)
{
	/* remove mutually exclusive command line options */
	StringT option(str);
	if (option == "-run")
	{
		RemoveCommandLineOption("-decomp");		
		RemoveCommandLineOption("-join");
	}
	else if (option == "-decomp")
	{
		RemoveCommandLineOption("-run");
		RemoveCommandLineOption("-join");
	}
	else if (option == "-join")
	{
		RemoveCommandLineOption("-decomp");		
		RemoveCommandLineOption("-run");
	}
	else if (option == "-join_io")
		RemoveCommandLineOption("-split_io");
	else if (option == "-split_io")
		RemoveCommandLineOption("-join_io");
	else if (option == "-decomp_method") {
	
		/* clear existing setting */
		int index;
		if (CommandLineOption("-decomp_method", index))
		{
			/* clear option */
			fCommandLineOptions.DeleteAt(index);

			/* clear number */
			const StringT& opt = fCommandLineOptions[index];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				fCommandLineOptions.DeleteAt(index);
		}
	}
	else if (option == "-dtd") {
		RunWriteDescription(XML_Attribute_FormatterT::DTD);
		return true;
	}
	else if (option == "-xsd") {
		RunWriteDescription(XML_Attribute_FormatterT::XSD);
		return true;
	}
	
	/* inherited */
	return ExecutionManagerT::AddCommandLineOption(option);
}

/* remove the command line option to the list */
bool FEExecutionManagerT::RemoveCommandLineOption(const char* str)
{
	StringT option(str);
	if (option == "-decomp")
	{
		int index;
		if (CommandLineOption("-decomp", index))
		{
			/* check for following number option */
			if (fCommandLineOptions.Length() > index+1) {
				const StringT& opt = fCommandLineOptions[index+1];
				if (strlen(opt) > 1 && isdigit(opt[1]))
					RemoveCommandLineOption(opt);
			}
			
			/* remove decomp */
			return ExecutionManagerT::RemoveCommandLineOption(option);
		}
		return false;
	}
	else if (option == "-join")
	{
		int index;
		if (CommandLineOption("-join", index))
		{
			/* check for following number option */
			if (fCommandLineOptions.Length() > index+1) {
				const StringT& opt = fCommandLineOptions[index+1];
				if (strlen(opt) > 1 && isdigit(opt[1]))
					RemoveCommandLineOption(opt);
			}
			
			/* remove decomp */
			return ExecutionManagerT::RemoveCommandLineOption(option);		
		}
		return false;
	} 
	else /* inherited */
		return ExecutionManagerT::RemoveCommandLineOption(option);
}

/**********************************************************************
 * Private
 **********************************************************************/

/* dump current DTD or XSD file */
void FEExecutionManagerT::RunWriteDescription(int doc_type) const
{
	try {
//TEMP - parameters currently needed to construct an FEManagerT
	StringT input;
	ofstreamT output;
	CommunicatorT comm;
	FEManagerT::TaskT task = FEManagerT::kParameters;
//TEMP

	/* collect parameters */
	cout << " collecting parameters..." << endl;
	ParameterTreeT tree;

	FEManagerT fe(input, output, comm, fCommandLineOptions, task);
	tree.BuildDescription(fe);

#ifdef BRIDGING_ELEMENT

	FEManagerT_bridging bridging(input, output, comm, fCommandLineOptions, task);
	tree.BuildDescription(bridging);

	MultiManagerT multi(input, output, comm, fCommandLineOptions, task);
	tree.BuildDescription(multi);

#ifdef __DEVELOPMENT__

	FEManagerT_THK thk(input, output, comm, fCommandLineOptions, task);
	tree.BuildDescription(thk);

	BridgingScaleManagerT bridging_multi(input, output, comm, fCommandLineOptions, task);
	tree.BuildDescription(bridging_multi);

#endif /* __DEVELOPMENT__ */
#endif /* BRIDGING_ELEMENT */

	/* write description */
	cout << " writing description..." << endl;
	StringT out_path("tahoe");
	XML_Attribute_FormatterT::DocTypeT doc_type_ = XML_Attribute_FormatterT::Undefined;
	if (doc_type == XML_Attribute_FormatterT::DTD) {
		doc_type_ = XML_Attribute_FormatterT::DTD;
		out_path.Append(".dtd");
	}
	else if (doc_type == XML_Attribute_FormatterT::XSD) {
		doc_type_ = XML_Attribute_FormatterT::XSD;
		out_path.Append(".xsd");
	}
	XML_Attribute_FormatterT attribute(doc_type_);

	ofstreamT out;
	out.open(out_path);
	attribute.InitDescriptionFile(out);
	
	const ArrayT<ParameterListT*>& branches = tree.Branches();
	for (int i = 0; i < branches.Length(); i++)
		attribute.WriteDescription(out, *(branches[i]));

	/* write statistics */
	cout << " " << attribute.ElementCount() << " elements" << '\n';
	cout << " " <<  attribute.AttributeCount() << " attributes" << '\n';
	cout << " " <<  attribute.LimitCount() << " limits" << '\n';

	attribute.CloseDescriptionFile(out);
	out.close();
	cout << " wrote \"" << out_path << '"' << endl;	
	}
	
	catch (ExceptionT::CodeT exc) {
		cout << "\n FEExecutionManagerT::RunDTD: caught exception: " 
		     << ExceptionT::ToString(exc) << endl;
	}
}

/* recursive dispatch */
void FEExecutionManagerT::JobOrBatch(ifstreamT& in, ostream& status)
{
	StringT ext;
	ext.Suffix(in.filename());
	if (ext == ".xml")
		RunJob(in, status);
	else /* inherited */
		ExecutionManagerT::JobOrBatch(in, status);
}

void FEExecutionManagerT::RunJob_analysis(const StringT& input_file, ostream& status) const
{
	const char* caller = "FEExecutionManagerT::RunJob_analysis";

	/* get rank and size */
	int rank = fComm.Rank();
	int size = fComm.Size();

	/* time markers */
	clock_t t0 = 0, t1 = 0, t2 = 0;
	
	/* start day/date info */
	time_t starttime;
	time(&starttime);

	/* output stream */
	StringT outfilename;
	outfilename.Root(input_file);
	if (size > 1) {outfilename.Append(".p", rank);}
	outfilename.Append(".out");
	ofstreamT out;
	out.open(outfilename);

#ifdef __MWERKS__
	if (!out.is_open()) {
		cout << "\n " << caller << ": could not open file: " << outfilename << endl;
		return;
	}
#endif

	int phase;
	int token = 1; /* for run time check sums */
	
	FEManagerT* tahoe = NULL;
	try
	{
		t0 = clock();
		phase = 0;

		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, true, true, fCommandLineOptions);

		/* write the validated list as formatted text */
		if (true) {
			DotLine_FormatterT pp_format;
			pp_format.SetTabWidth(4);
			pp_format.InitParameterFile(out);
			pp_format.WriteParameterList(out, valid_list);
			pp_format.CloseParameterFile(out);
			out << endl;
		}

		/* construction */
		tahoe = FEManagerT::New(valid_list.Name(), input_file, out, fComm, fCommandLineOptions, FEManagerT::kRun);
		
		tahoe->TakeParameterList(valid_list);

		if (fComm.Sum(token) != size && size > 1)
		{
			/* gather tokens to rank 0 */
			if (rank == 0)
			{
				iArrayT tokens(size);
				fComm.Gather(token, tokens);
			
				status << "\n The following processes exit on exception during construction:\n";
				for (int i = 0; i < tokens.Length(); i++)
					if (tokens[i] != 1)
						status << setw(kIntWidth) << i << '\n';
				status.flush();
			}
			else fComm.Gather(token, 0);
			
			ExceptionT::GeneralFail(caller);
		}

		IOManager_mpi* IOMan = NULL;
		if (fComm.Size() > 1 && valid_list.Name() == "tahoe")
		{
			/* partition information */
			const PartitionT* partition = tahoe->Partition();
		
			/* external IO */
			token = 1;
			if (partition && partition->DecompType() == PartitionT::kGraph && !CommandLineOption("-split_io"))
			{
				try {

					/* geometry file information */
					const StringT& model_file = valid_list.GetParameter("geometry_file");
					IOBaseT::FileTypeT format = IOBaseT::int_to_FileTypeT(valid_list.GetParameter("geometry_format"));
					
					/* set-up local IO */
					IOMan = new IOManager_mpi(input_file, fComm, *(tahoe->OutputManager()), *partition, model_file, format);
					if (!IOMan) throw ExceptionT::kOutOfMemory;
				
					/* set external IO */
					tahoe->SetExternalIO(IOMan);
				}
			
				catch (ExceptionT::CodeT code)
				{
					token = 0;
					status << "\n \"" << input_file << "\" exit on exception " << code 
				    	   << " setting the external IO" << endl;
				}
			}
		}
		if (fComm.Sum(token) != size && size > 1)
		{
			/* gather tokens to rank 0 */
			if (rank == 0)
			{
				iArrayT tokens(size);
				fComm.Gather(token, tokens);			
				status << "\n The following processes exit on exception during construction:\n";
				for (int i = 0; i < tokens.Length(); i++)
					if (tokens[i] != 1)
						status << setw(kIntWidth) << i << '\n';
				status.flush();
			}
			else 
				fComm.Gather(token, 0);
			
			ExceptionT::GeneralFail(caller);
		}
		
		
#if defined(__MWERKS__) && __option(profile)
		/* start recording profiler information */
		ProfilerSetStatus(1);
#endif
		
		/* solution */
		phase = 1;
		t1 = clock();
		try 
		{
			tahoe->Solve();
		}
		catch(ExceptionT::CodeT code)
		{
			status << "\n \"" << input_file << "\" exit on exception " << code << " during the\n";
			status << " solution phase. See \"" << outfilename << "\" for a list";
			status << " of the codes.\n";
			token = 0;
		}
		t2 = clock();

#if defined(__MWERKS__) && __option(profile)
		/* stop recording profiler information */
		ProfilerSetStatus(0);
#endif
		if (valid_list.Name() == "tahoe")
		{
			/* free external IO */
			delete IOMan;
		}
		if (fComm.Sum(token) != size && size > 1)
		{
			/* gather tokens to rank 0 */
			if (rank == 0)
			{
				iArrayT tokens(size);
				fComm.Gather(token, tokens);

				status << "\n The following processes exit on exception during solution:\n";
				for (int i = 0; i < tokens.Length(); i++)
					if (tokens[i] != 1)
						status << setw(kIntWidth) << i << '\n';
				status.flush();
			}
			else fComm.Gather(token, 0);
		
			ExceptionT::GeneralFail(caller);
		}
		/* clean up */
		delete tahoe;
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		/* clean up */
		delete tahoe;
	
		status << "\n \"" << input_file << "\" exit on exception during the\n";
		if (phase == 0)
		{
			status << " construction phase. Check the input file for errors." << endl;
		}
		else
		{
			status << " solution phase. See \"" << outfilename << "\" for a list";
			status << " of the codes.\n";
		}
		
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		

		out << endl;
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	
	out << "\n   Start time: " << ctime(&starttime);
	out <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "     Solution: "   << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status << "    Stop time: " << ctime(&stoptime);
	out   << "    Stop time: " << ctime(&stoptime);

	status << "\n End Execution\n" << endl;
	out    << "\n End Execution\n" << endl;
}	

/* generate decomposition files */
void FEExecutionManagerT::RunDecomp_serial(const StringT& input_file, ostream& status,CommunicatorT& comm, int size) const
{
	/* look for size */
	int index;
	if (CommandLineOption("-decomp", index) && fCommandLineOptions.Length() > index+1) {
		const char* opt = fCommandLineOptions[index+1];
		if (strlen(opt) > 1 && isdigit(opt[1]))
			size = atoi(opt+1); /* opt[0] = '-' */
	}

	/* look for method */
	int method = -1;
	if (CommandLineOption("-decomp_method", index))
	{	
		if (fCommandLineOptions.Length() > index+1)
		{
			const char* opt = fCommandLineOptions[index+1];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				method = atoi(opt+1); /* opt[0] = '-' */
		}
	}
	
	/* prompt for size */
	if (size == -1) 
	{
		/* prompt for decomp size */
		int count = 0;
		while (count == 0 || (count < 5 && size < 2))
		{
			count++;

			/* number of partitions */					
			cout << "\n Enter number of partitions > 1 (0 to quit): ";
#if (defined __SGI__ && defined __TAHOE_MPI__)
			cout << '\n';
#endif
			cin >> size;
		
			/* clear to end of line */
			fstreamT::ClearLine(cin);

			if (size == 0) break;
		}
	}
	if (size < 2) return;

	/* prompt for method */
	if (method == -1)
	{
		cout << "\n Select partitioning method:\n"
		     << '\t' << PartitionT::kGraph   << ": graph\n"
		     << '\t' << PartitionT::kAtom    << ": atom\n"
		     << '\t' << PartitionT::kSpatial << ": spatial\n";
		cout << "\n method: "; 
#if (defined __SGI__ && defined __TAHOE_MPI__)
		cout << '\n';
#endif					
		cin >> method;
	
		/* clear to end of line */
		fstreamT::ClearLine(cin);
	}

	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();

		/* path to parameters file */
		StringT path;
		path.FilePath(input_file);
		
		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, false, false, fCommandLineOptions);
		
		/* extract model file and file format */
		int i_format = valid_list.GetParameter("geometry_format");
		IOBaseT::FileTypeT format = IOBaseT::int_to_FileTypeT(i_format);
		StringT model_file = valid_list.GetParameter("geometry_file");

		/* name translation */
		model_file.ToNativePathName();      
		model_file.Prepend(path);
		
		/* set output map and and generate decomposition */
		DecomposeT decompose;
		decompose.CheckDecompose(input_file, size, method, comm, model_file, format, fCommandLineOptions);

		t1 = clock();
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
}

/* join parallel results files */
void FEExecutionManagerT::RunJoin_serial(const StringT& input_file, ostream& status, int size) const
{
	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();

		/* path to parameters file */
		StringT path;
		path.FilePath(input_file);

		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, false, false, fCommandLineOptions);
		
		/* model file parameters */
		int i_format = valid_list.GetParameter("geometry_format");
		IOBaseT::FileTypeT model_format = IOBaseT::int_to_FileTypeT(i_format);
		i_format = valid_list.GetParameter("output_format");
		IOBaseT::FileTypeT results_format = IOBaseT::int_to_FileTypeT(i_format);
		StringT model_file = valid_list.GetParameter("geometry_file");
		if (results_format == IOBaseT::kTahoe ||
		    results_format == IOBaseT::kTahoeII)
			results_format = IOBaseT::kTahoeResults;

		/* name translation */
		model_file.ToNativePathName();      
		model_file.Prepend(path);

		int index;
		if (CommandLineOption("-join", index) && fCommandLineOptions.Length() > index+1) {
			const char* opt = fCommandLineOptions[index+1];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				size = atoi(opt+1);
		}

		/* prompt for decomp size */
		if (size == -1)
		{
			cout << "\n Enter number of partitions > 1 (0 to quit): ";
#if (defined __SGI__ && defined __TAHOE_MPI__)
			cout << '\n';
#endif					
			cin >> size;

			/* clear to end of line */
			fstreamT::ClearLine(cin);
		}

		if (size < 2) return;
		
		/* construct output formatter */
		StringT program = "tahoe";
		StringT version = FEManagerT::Version();
		StringT title;
		const ParameterT* title_param = valid_list.Parameter("title");
		if (title_param)
			title = *title_param;
		StringT input = input_file;
		OutputBaseT* output = IOBaseT::NewOutput(program, version, title, input, 
			results_format, cout);

		/* construct joiner */
		JoinOutputT output_joiner(input_file, model_file, model_format, 
			results_format, output, size);

		/* join files */
		output_joiner.Join();
		
		/* free output formatter */
		delete output;

		t1 = clock();
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "         Join: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
}
