/* $Id: main.cpp,v 1.11 2002-08-15 16:49:57 sawimme Exp $ */
/* created: paklein (05/22/1996) */
#include <iostream.h>
#include <fstream.h>

#include "Environment.h"
#include "ExceptionCodes.h"

#ifdef __MWERKS__
#if __option(profile)
#include <profiler.h>
#endif
#ifdef macintosh
extern "C" int ccommand(char ***arg);
#endif
#endif

/* MP environment */
#include "CommunicatorT.h"

#include "FEExecutionManagerT.h"
#include "StringT.h"

using namespace Tahoe;

static void StartUp(int* argc, char*** argv, CommunicatorT& comm);
static void ShutDown(CommunicatorT& comm);

/* redirect of cout for parallel execution */
ofstream console;
#if defined (__DEC__) || defined (__SUN__)
streambuf* cout_buff = NULL,*cerr_buff = NULL;
#endif

/* f2c library global variables */
int xargc;
char **xargv;

void main(int argc, char* argv[])
{
	/* f2c library global variables */
	xargc = argc;
	xargv = argv;
	
	/* MP */
	CommunicatorT::SetArgv(&argc, &argv);
	CommunicatorT comm;

	StartUp(&argc, &argv, comm);

	FEExecutionManagerT FEExec(argc, argv, '%', '@', comm);
	FEExec.Run();		

	ShutDown(comm);
}

static void StartUp(int* argc, char*** argv, CommunicatorT& comm)
{
#if !defined(_MACOS_)
#if defined (__DEC__) || defined (__SUN__)
	/* redirect cout and cerr */
	if (comm.Rank() > 0)
	{
		StringT console_file("console");
		console_file.Append(comm.Rank());
		console.open(console_file, ios::app);

		/* keep buffers from cout and cerr */
		cout_buff = cout.rdbuf();
		cerr_buff = cerr.rdbuf();

		/* redirect */
		cout.rdbuf(console.rdbuf());
		cerr.rdbuf(console.rdbuf());
	}
#else
	/* redirect cout and cerr */
	if (comm.Rank() > 0)
	{
		StringT console_file("console");
		console_file.Append(comm.Rank());
		console.open(console_file, ios::app);
		cout = console;
		cerr = console;
	}
#endif /* __DEC__ */
#else /* __MACOS__ */
#pragma unused(comm)
#endif /* __MACOS__ */

	/* output build date and time */
	cout << "\n build: " __TIME__ ", " << __DATE__ << '\n';

#if defined(__MWERKS__) && defined (macintosh)
	/* get command-line arguments */
	*argc = ccommand(argv);
#endif /* Carbon */

	if (CommunicatorT::ActiveMP())
		cout << "********************* MPI Version *********************" << '\n';

#if __option (extended_errorcheck)
	cout << "\n Extended error checking is ON\n";
	cout << " Turn off \"extended error checking\" to remove.\n";
#else
	cout << "\n Extended error checking is OFF\n";
	cout << " Turn on \"extended error checking\" to add.\n";
#endif

#if __option (profile)
	ProfilerInit(collectDetailed,bestTimeBase,100,20);
	ProfilerSetStatus(0);
	cout << "\n P r o f i l i n g   i s   a c t i v e\n" << endl;
#endif

#ifdef __NO_RTTI__
	cout << "\n RTTI not supported. Not all options will be\n";
	cout << " available.\n" << endl;
#endif
}

static void ShutDown(CommunicatorT& comm)
{
	cout << "\nExit.\n" << endl;

#if __option (profile)		
	ProfilerDump("\ptahoe.prof");
	ProfilerTerm();
#endif

	if (comm.Rank() > 0)
	{
#ifdef __DEC__
		/* restore cout and cerr */
		cout.rdbuf(cout_buff);
		cerr.rdbuf(cerr_buff);
#endif
		/* close console file */
		console.close();
	}
}
