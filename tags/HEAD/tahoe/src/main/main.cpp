/* $Id: main.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

#include <iostream.h>
#include <fstream.h>

#include "Environment.h"
#include "ExceptionCodes.h"

#if defined(__MWERKS__) && __option(profile)
#include <profiler.h>
#endif

#ifdef __MPI__
#include "mpi.h"
#endif

#include "FEExecutionManagerT.h"
#include "StringT.h"

static void StartUp(int* argc, char*** argv);
static void ShutDown(void);

/* redirect of cout for parallel execution */
ofstream console;

/* f2c library global variables */
int xargc;
char **xargv;

void main(int argc, char* argv[])
{
	/* f2c library global variables */
	xargc = argc;
	xargv = argv;

	StartUp(&argc, &argv);

	FEExecutionManagerT FEExec(argc, argv, '%', '@');
	FEExec.Run();		

	ShutDown();
}

static void StartUp(int* argc, char*** argv)
{

#ifdef __MPI__
#ifdef __MWERKS__
	cout << " Initializing MPI..." << endl;
#endif

	/* initialize MPI */
	if (MPI_Init(argc, argv) != MPI_SUCCESS)
	{
		cout << "\n MPI_Init: error" << endl;
		throw eMPIFail;
	}

	int rank;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
		throw eMPIFail;
	
#ifndef _MACOS_
	/* redirect cout and cerr */
	if (rank > 0)
	{
		StringT console_file("console");
		console_file.Append(rank);
		console.open(console_file, ios::app);
		cout = console;
		cerr = console;
	}
#endif /* __MACOS__ */
#else
#pragma unused(argc)
#pragma unused(argv)
#endif /* __MPI__ */

#ifndef _MACOS_
	cout << "********************* MPI Version *********************" << '\n';
#endif /* _MACOS_ */

	/* output build date and time */
	cout << "\n build: " __TIME__ ", " << __DATE__ << '\n';

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

static void ShutDown(void)
{
	cout << "\nExit.\n" << endl;

#if __option (profile)		
	ProfilerDump("\ptahoe.prof");
	ProfilerTerm();
#endif

#ifdef __MPI__
	/* get rank */
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* end MPI */
	if (MPI_Finalize() != MPI_SUCCESS)
		cout << "\n MPI_Finalize: error" << endl;

	if (rank > 0) console.close();
#endif /* __MPI__ */
}
