/* $Id: ExecutionManagerT.cpp,v 1.1.1.1 2001-01-29 08:20:21 paklein Exp $ */
/* created: paklein (08/27/1997)                                          */
/* Manages input file driven jobs.                                        */
/* MUST overload private:RunJob().                                        */

#include "ExecutionManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>

#ifdef __MPI__
#include "mpi.h"
#endif

#include "fstreamT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "StringT.h"

/* maximum batch file recursion depth */
const int kMaxRecursionDepth = 10;

/* Constructor */
ExecutionManagerT::ExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
	int jobcharputback):
	fJobChar(job_char),
	fBatchChar(batch_char),
	fJobCharPutBack(jobcharputback),
	fRecursionDepth(0)
{
	if (fJobCharPutBack != 1 && fJobCharPutBack != 0) throw eBadInputValue;

	/* store command line arguments */
	fCommandLineOptions.Allocate(argc);
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
		fCommandLineOptions[i] = argv[i];

//TEMP - no command line arguments for Mac
#ifdef _MACOS_	
	bool decomp_test = false;
	if (decomp_test)
	{
		cout << "\n ExecutionManagerT::ExecutionManagerT: command-line test" << endl;
		fCommandLineOptions.Allocate(2);
		fCommandLineOptions[0] = "tahoe";
		fCommandLineOptions[1] = "-decomp";
	}
#endif

	//TEMP
	if (0 && fCommandLineOptions.Length() > 1)
	{
		cout << "\n command line arguments:\n";
		for (int j = 0; j < fCommandLineOptions.Length(); j++)
			cout << setw(kIntWidth) << j << ":" << fCommandLineOptions[j] << '\n';
		cout.flush();
	}

	/* format standard output */
	ofstreamT::format_stream(cout);	
}

/* Destructor */
ExecutionManagerT::~ExecutionManagerT(void) { }

/* Prompt input files until "quit" */
void ExecutionManagerT::Run(void)
{
	int size = 1;
#ifdef __MPI__
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;
#endif

	/* dispatch */
	if (size > 1)
		Run_parallel();
	else
		Run_serial();
}

/**********************************************************************
* Protected
**********************************************************************/

void ExecutionManagerT::Run_serial(void)
{
	/* file name passed on command line */
	int index;
	if (CommandLineOption("-f", index))
	{
		/* path name */
		StringT& file = fCommandLineOptions[index+1];
		file.ToNativePathName();
	
		/* open stream */
		ifstreamT input(file);
		if (!input.is_open())
		{
			cout << "\n ExecutionManagerT::Run_serial: unable to open file: \""
			     << file  << '\"' << endl;
			throw eBadInputValue;
		}
		
		/* dispatch */
		JobOrBatch(input, cout);
	}
	else
	{
		StringT lastfilename;
		while (1)
		{
			/* prompt for input filename and open stream */
			ifstreamT input;		
			if (!input.open("Enter input file name", "quit", lastfilename)) break;
			
			/* keep last file name */
			lastfilename = input.filename();
	
			/* Recursive dispatch */
			JobOrBatch(input, cout);
		}
	}
}

void ExecutionManagerT::Run_parallel(void)
#ifdef __MPI__
{
	/* get rank */
	int rank;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) throw eMPIFail;

	/* file name passed on command line */
	int index;
	if (CommandLineOption("-f", index))
	{
		/* path name */
		StringT& file = fCommandLineOptions[index+1];
		file.ToNativePathName();

		/* open stream */
		ifstreamT input(file);
		if (!input.is_open())
		{
			cout << "\n ExecutionManagerT::Run_parallel: unable to open file: \""
			     << file  << '\"' << endl;
			throw eBadInputValue;
		}
		
		/* dispatch */
		JobOrBatch(input, cout);
	}
	else
	{

		StringT lastfilename;
		while (1)
		{
			/* prompt for input filename and open stream */
			ifstreamT input;
			StringT file(255);
			if (rank == 0)
			{
				if (input.open("Enter input file name", "quit", lastfilename))
					file.CopyIn(input.filename());
				else
					file = "quit";
			}
	
			/* broadcast file name */
			if (MPI_Bcast(file.Pointer(), 255, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
				throw eMPIFail;		
			
			/* open stream */
			if (file == "quit")
				break;
			else if (rank > 0)
				input.open(file);
		
			/* keep last file name */
			lastfilename = input.filename();
	
			/* recursive dispatch */
			JobOrBatch(input, cout);
		}
	}
}
#else /* __MPI__ */
{
	cout << "\n ExecutionManagerT::Run_parallel: should not be here without MPI" << endl;
	throw eGeneralFail;
}
#endif /* __ MPI__ */

/* returns true if the option was passed on the command line */
bool ExecutionManagerT::CommandLineOption(const char* str, int& index) const
{
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
		if (fCommandLineOptions[i] == str)
		{
			index = i;
			return true;
		}

	/* dummy */
	index = 0;
	return false;
}

void ExecutionManagerT::AddCommandLineOption(const char* str)
{
	/* only if not present */
	int index;
	if (!CommandLineOption(str, index))
	{
		int num_options = fCommandLineOptions.Length();
	
		ArrayT<StringT> temp(num_options + 1);
		for (int i = 0; i < num_options; i++)
			temp[i] = fCommandLineOptions[i];
	
		/* add new */
		temp[num_options] = str;
		
		/* exhange data */
		fCommandLineOptions.Swap(temp);
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* Recursive dispatch */
void ExecutionManagerT::JobOrBatch(ifstreamT& in, ostream& status)
{
	/* get first char */
	char filetypechar;
	in >> filetypechar;
	
	if (filetypechar != fJobChar && filetypechar != fBatchChar)
	{
		status << "\n ExecutionManagerT::JobOrBatch: invalid filetype character: ";
		status << filetypechar << '\n';
		return;
	}

	/* check recursion depth */
	if (++fRecursionDepth > kMaxRecursionDepth) throw eGeneralFail;
	
	/* JOB file */
	if (filetypechar == fJobChar)
	{
		/* putback first character */
		if (fJobCharPutBack) in.putback(filetypechar);
		
		/* derived */
		RunJob(in, status);
	}
	/* BATCH file */
	else
	{
		/* process batch file */
		RunBatch(in, status);
	}
	
	/* reduce depth on exit */
	fRecursionDepth--;
}
	
/* Batch file processing */
void ExecutionManagerT::RunBatch(ifstreamT& in, ostream& status)
{
#ifdef __MPI__
	int size;
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;
	if (size > 1)
	{
		cout << "\n ExecutionManagerT::RunBatch: not ready for parallel execution" << endl;
		throw eGeneralFail;
	}
#endif /* __MPI__ */

	/* mark status */
	status << "\n Processing batch file: " << in.filename() << '\n';

	/* clear whitespace */
	//in.next_char();

	/* open status stream */
	StringT statusfilename;
	ofstreamT stat;	
	stat.open(statusfilename.DefaultName(in.filename(),".bat",".stat", -1));
	
	/* start day/date info */
	time_t starttime;
	time(&starttime);

	/* get 1st entry */
	StringT nextinfilename;
	in >> nextinfilename;
	
	/* repeat to end of file */
	while (in.good())
	{
		/* open new input stream */
		nextinfilename.ToNativePathName();
		ifstreamT nextin(nextinfilename);
	
		/* process if valid */
		if (nextin.is_open())
			JobOrBatch(nextin, stat);
		else
			stat << " File not found: " << nextinfilename << '\n';
			
		/* get next entry */
		in >> nextinfilename;
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);
	stat << "\n Batch start time  : " << ctime(&starttime);
	stat <<   " Batch stop time   : " << ctime(&stoptime);
}
