/* $Id: main.cpp,v 1.22.26.1 2004-07-13 16:42:41 paklein Exp $ */
/* created: paklein (05/22/1996) */
#include <iostream.h>
#include <fstream.h>

#include "Environment.h"
#include "ExceptionT.h"

#ifdef __MWERKS__
#if __option(profile)
#include <Profiler.h>
#endif
#ifdef macintosh
extern "C" int ccommand(char ***arg);
#endif
#endif

#if defined(__MWERKS__) && defined(__MACH__)
#include <unistd.h> //TEMP - needed for chdir
#endif

/* MP environment */
#include "CommunicatorT.h"

#include "FEExecutionManagerT.h"
#include "StringT.h"

using namespace Tahoe;

static void StartUp(int* argc, char*** argv, CommunicatorT& comm);
static void ShutDown(CommunicatorT& comm);
static void DumpLicense(void);

/* redirect of cout for parallel execution */
ofstream console;
#if defined (__DEC__) || defined (__SUN__) || defined(__GCC_3__) || defined(__INTEL_CC__)
streambuf* cout_buff = NULL,*cerr_buff = NULL;
#endif

/* f2c library global variables */
int xargc;
char **xargv;

int main(int argc, char* argv[])
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
	return 0;
}

static void StartUp(int* argc, char*** argv, CommunicatorT& comm)
{
#if !defined(_MACOS_) && !defined(__INTEL__)
#if defined (__DEC__) || defined (__SUN__) || defined(__GCC_3__) || defined(__INTEL_CC__)
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
#endif /* __DEC__ || __SUN__ || __GCC_3__ || __INTEL_CC__ */
#else /* __MACOS__ && __INTEL__ */
#pragma unused(comm)
#endif /* __MACOS__ && __INTEL__ */

	/* write license to screen */
	DumpLicense();	

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
	pascal OSErr err = ProfilerInit(collectDetailed, bestTimeBase, 1000, 25);
	if (err == memFullErr)
		cout << "\n ProfilerInit: memory full" << endl;
	else if	(err == paramErr)
		cout << "\n ProfilerInit: parameter error" << endl;
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

void DumpLicense(void)
{
	const char version[] = "Tahoe 2.0";
	cout << "\n " << version << "\n\n"
         << " Copyright 2003, Sandia Corporation.\n" 
	     << " All rights reserved.\n\n"
	     << " Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive\n"
	     << " license for use of this work by or on behalf of the U.S. Government. Export\n"
	     << " of this program may require a license from the United States Government."
	     << endl;
}
