/* $Id: ExecutionManagerT.cpp,v 1.17.2.1 2004-09-15 02:14:15 d-farrell2 Exp $ */
/* created: paklein (08/27/1997) */
#include "ExecutionManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
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
		int count = 0;
		while (count++ < 10)
		{
			/* broadcast file name or command-line option */
			StringT line(255);
			if (rank == 0)
				Prompt("Enter input file path or option (\"quit\" to exit)", lastfilename, line);
			fComm.Broadcast(0, line);
		
			/* command line option */
			if (line[0] == '-') {

				/* reset count */
				count = 0;

				/* add command line option */
				if (AddCommandLineOption(line))
					cout << " added command line option: \"" << line << '\"' << endl;
			}
			else if (strncmp(line, "quit", 4) == 0)
				break;
			else /* try to open file */
			{
				line.ToNativePathName();
				ifstreamT input('#', line);
				if (input.is_open()) {

					/* reset count */
					count = 0;

					/* keep last file name */
					lastfilename = input.filename();
	
					/* Recursive dispatch */
					JobOrBatch(input, cout);
				}
				else
					cout << "\nError: filename: \"" << line << "\" not found\n";
			}
		}
		
		/* exit message */
		if (count >= 10) cout << "\nNo valid input after " << count << " attempts." <<  endl;
	}
}

// returns true if the option was passed on the command line, moved to public DEF 3 Aug 04
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


/**********************************************************************
 * Protected
 **********************************************************************/

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

void ExecutionManagerT::Prompt(const char* prompt, const char* default_input, StringT& line) const
{
	cout << '\n' << prompt;

	/* default */
	if (default_input != NULL && strlen(default_input) > 0)
	{
		cout << "\nEnter <RETURN> for \"" << default_input << "\": ";
#ifdef __SGI__
		cout.flush();
#endif
			
		/* new filename */
		char test = cin.peek();
		if (test != '\n')
		{
			/* take first word */
			cin >> line;
		}
		else
		{
			/* copy default */
			line = default_input;
		}				
	}
	else
	{
		cout << ": ";
#ifdef __SGI__
		cout.flush();
#endif					
		cin >> line;
	}
		
	/* clear to end of line */
	fstreamT::ClearLine(cin);
}
