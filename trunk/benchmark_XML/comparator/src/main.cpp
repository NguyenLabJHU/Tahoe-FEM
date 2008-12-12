/* $Id: main.cpp,v 1.6 2008-12-12 17:50:43 lxmota Exp $ */
/* created: paklein (05/22/1996) */
#include <iostream.h>
#include <fstream.h>

#include "Environment.h"

#ifdef __MWERKS__
#if __option(profile)
#include <profiler.h>
#endif
#ifdef macintosh
extern "C" int ccommand(char ***arg);
#endif
#endif

#include "ComparatorT.h"

using namespace Tahoe;

static void StartUp(int* argc, char*** argv);
static void ShutDown(void);

int main(int argc, char* argv[])
{
	StartUp(&argc, &argv);

	ComparatorT compare(argc, argv, '%', '@');
	compare.Run();

	ShutDown();
	return 0;
}

static void StartUp(int* argc, char*** argv)
{
#pragma unused(argc)
#pragma unused(argv)

#if defined(__MWERKS__) && defined (macintosh)
	/* get command-line arguments */
	*argc = ccommand(argv);
#endif

#if __option (extended_errorcheck)
	cout << "\n Extended error checking is ON\n";
	cout << " Turn off \"extended error checking\" to remove.\n";
#else
	cout << "\n Extended error checking is OFF\n";
	cout << " Turn on \"extended error checking\" to add.\n";
#endif

#if __option (profile)
	ProfilerInit(collectDetailed, bestTimeBase, 100, 20);
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
}
