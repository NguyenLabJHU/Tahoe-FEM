/* $Id: FileCrawlerT.cpp,v 1.7 2004-06-17 06:38:29 paklein Exp $ */

#include "FileCrawlerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "StringT.h"

/* maximum batch file recursion depth */

using namespace Tahoe;

const int kMaxRecursionDepth = 10;

/* Constructor */
FileCrawlerT::FileCrawlerT(int argc, char* argv[], char job_char, char batch_char):
	fJobChar(job_char),
	fBatchChar(batch_char),
	fRecursionDepth(0)
{
	/* store command line arguments */
	fCommandLineOptions.Dimension(argc);
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
		fCommandLineOptions[i] = argv[i];

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
FileCrawlerT::~FileCrawlerT(void) { }

/* Prompt input files until "quit" */
void FileCrawlerT::Run(void)
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
			cout << "\n FileCrawlerT::Run: unable to open file: \""
			     << file  << '\"' << endl;
			throw ExceptionT::kBadInputValue;
		}
		
		/* dispatch */
		try { JobOrBatch(input, cout); }
		catch (ExceptionT::CodeT error)
		{
			cout << "\n FileCrawlerT::Run: file \"" << file 
			     << "\" quit on exception: " << error << endl;
		}
	}
	else
	{
		StringT lastfilename;
		while (1)
		{
			/* prompt for input filename and open stream */
			ifstreamT input('#');
			if (!input.open("Enter input file name", "quit", lastfilename)) break;
			
			/* keep last file name */
			lastfilename = input.filename();
	
			/* Recursive dispatch */
			try { JobOrBatch(input, cout); }
			catch (ExceptionT::CodeT error)
			{
				cout << "\n FileCrawlerT::Run: file \"" << input.filename() 
				     << "\" quit on exception: " << error << endl;
			}
		}
	}
}

/**********************************************************************
* Protected
**********************************************************************/

/* Batch file processing */
void FileCrawlerT::RunBatch(ifstreamT& in, ostream& status)
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
		/* command in batch file */
		if (nextinfilename[0] == '-')
		{
			BatchFileCommand(nextinfilename, in);
		}
		else
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

/* returns true if the option was passed on the command line */
bool FileCrawlerT::CommandLineOption(const char* str, int& index) const
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

void FileCrawlerT::AddCommandLineOption(const char* str)
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
void FileCrawlerT::JobOrBatch(ifstreamT& in, ostream& status)
{
	/* get first char */
	char filetypechar;
	in >> filetypechar;
	
	if (filetypechar != fJobChar && filetypechar != fBatchChar)
	{
		status << "\n FileCrawlerT::JobOrBatch: invalid filetype character: ";
		status << filetypechar << '\n';
		return;
	}

	/* check recursion depth */
	if (++fRecursionDepth > kMaxRecursionDepth) throw ExceptionT::kGeneralFail;
	
	/* JOB file */
	if (filetypechar == fJobChar)
	{
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
