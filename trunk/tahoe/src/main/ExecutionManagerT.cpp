/* $Id: ExecutionManagerT.cpp,v 1.14 2003-01-27 07:00:27 paklein Exp $ */
/* created: paklein (08/27/1997) */
#include "ExecutionManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>

#include "fstreamT.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "StringT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* maximum batch file recursion depth */
const int kMaxRecursionDepth = 10;

/* Constructor */
ExecutionManagerT::ExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
	CommunicatorT& comm,
	int jobcharputback):
	fJobChar(job_char),
	fBatchChar(batch_char),
	fComm(comm),
	fJobCharPutBack(jobcharputback),
	fRecursionDepth(0)
{
	if (fJobCharPutBack != 1 && fJobCharPutBack != 0) throw ExceptionT::kBadInputValue;

	/* store command line arguments */
	fCommandLineOptions.Dimension(argc);
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
		fCommandLineOptions[i] = argv[i];

	/* format standard output */
	ofstreamT::format_stream(cout);	
}

/* Destructor */
ExecutionManagerT::~ExecutionManagerT(void) { }

/* Prompt input files until "quit" */
void ExecutionManagerT::Run(void)
{
	/* dispatch */
	if (fComm.Size() > 1)
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
		ifstreamT input('#', file);
		if (!input.is_open())
		{
			cout << "\n ExecutionManagerT::Run_serial: unable to open file: \""
			     << file  << '\"' << endl;
			throw ExceptionT::kBadInputValue;
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
			ifstreamT input('#');
			if (!OpenWithPrompt("Enter input file name or command line option", "quit", lastfilename, input)) break;
			
			/* keep last file name */
			lastfilename = input.filename();
	
			/* Recursive dispatch */
			JobOrBatch(input, cout);
		}
	}
}

void ExecutionManagerT::Run_parallel(void)
{
	/* get rank */
	int rank = fComm.Rank();

	/* file name passed on command line */
	int index;
	if (CommandLineOption("-f", index))
	{
		/* path name */
		StringT& file = fCommandLineOptions[index+1];
		file.ToNativePathName();

		/* open stream */
		ifstreamT input('#', file);
		if (!input.is_open())
		{
			cout << "\n ExecutionManagerT::Run_parallel: unable to open file: \""
			     << file  << '\"' << endl;
			throw ExceptionT::kBadInputValue;
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
			ifstreamT input('#');
			StringT file(255);
			if (rank == 0)
			{
				if (input.open("Enter input file name", "quit", lastfilename))
					file.CopyIn(input.filename());
				else
					file = "quit";
			}
	
			/* broadcast file name */
			fComm.Broadcast(0, file);
			
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

bool ExecutionManagerT::AddCommandLineOption(const char* str)
{
	/* only if not present */
	int index;
	if (!CommandLineOption(str, index))
	{
		/* append */
		fCommandLineOptions.Append(str);
		return true;
	}
	else return false;
}

bool ExecutionManagerT::RemoveCommandLineOption(const char* str)
{
	/* only if not present */
	int index;
	if (CommandLineOption(str, index))
	{
		/* append */
		fCommandLineOptions.DeleteAt(index);
		return true;
	}
	else return false;
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
	if (++fRecursionDepth > kMaxRecursionDepth) throw ExceptionT::kGeneralFail;
	
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
	/* mark status */
	status << "\n Processing batch file: " << in.filename() << '\n';
	
	/* start day/date info */
	time_t starttime;
	time(&starttime);

	/* get 1st entry */
	StringT nextinfilename;
	in >> nextinfilename;
	
	/* repeat to end of file */
	while (in.good())
	{
		/* adjusting execution options */
		if (nextinfilename[0] == '-')
			AddCommandLineOption(nextinfilename);
		else /* execute regular file */
		{	
			/* file path format */
			nextinfilename.ToNativePathName();

			/* path to source file */
			StringT path;
			path.FilePath(in.filename());
	
			/* open new input stream */
			nextinfilename.Prepend(path);
			ifstreamT nextin('#', nextinfilename);
	
			/* process if valid */
			if (nextin.is_open())
				JobOrBatch(nextin, cout);
			else
				cout << " File not found: " << nextinfilename << '\n';
		}
			
		/* get next entry */
		in >> nextinfilename;
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);
	cout << "\n Batch start time  : " << ctime(&starttime);
	cout <<   " Batch stop time   : " << ctime(&stoptime);
}

/* open stream with prompt - return 1 if successful */
int ExecutionManagerT::OpenWithPrompt(const char* prompt, const char* skipname,
	const char* defaultname, ifstreamT& in)
{
	StringT newfilename;
	int maxtry = 20;
	int count  = 0;
	while (1)
	{
		cout << '\n' << prompt << "\n(\"" << skipname << "\" to exit";

		/* default */
		if (defaultname != NULL && strlen(defaultname) > 0)
		{
			cout << ", <RETURN> for \"" << defaultname << "\"): ";
#ifdef __SGI__
			cout.flush();
#endif
			
			/* new filename */
			char test = cin.peek();
			if (test != '\n')
			{
				/* take first word */
				cin >> newfilename;
			}
			else
			{
				/* copy default */
				newfilename = defaultname;
			}				
		}
		else
		{
			cout << "): ";
#ifdef __SGI__
			cout.flush();
#endif					
			cin >> newfilename;
		}
		
		/* clear to end of line */
		fstreamT::ClearLine(cin);

		/* check exit */
		if (strcmp(newfilename, skipname) == 0)
			return 0;
		/* check program option */
		else if (newfilename[0] == '-')
		{
			if (AddCommandLineOption(newfilename))
				cout << " added command line option: \"" << newfilename << '\"' << endl;
		}	
		/* attempt open file */
		else
		{
			/* convert to native path */
			newfilename.ToNativePathName();
		
			/* attempt open */
			in.open(newfilename);
			
			/* check file open */
			if (in.is_open())
				return 1;	
			else
			{
				cout << "\nError: filename: " << newfilename << " not found\n";
				in.clear();
			}
			
			/* could not find file */
			if (++count == maxtry)
			{
				cout << "\n ExecutionManagerT::OpenInputStream: could not find file after ";
				cout << maxtry << " iterations" << endl;
				throw ExceptionT::kGeneralFail;
			}
		}
	}	
}

