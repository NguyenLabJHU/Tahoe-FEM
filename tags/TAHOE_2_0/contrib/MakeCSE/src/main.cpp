// $Id: main.cpp,v 1.10 2003-09-05 19:48:55 paklein Exp $
// created: 6 Oct 1999 by S. A. Wimmer
// program reads input file, runs MakeCSE, writes output file

#include "ExceptionT.h"
#include "MakeCSE_ExecutionT.h"
#include "sArrayT.h"

#if defined(__MWERKS__) && defined(__MACH__)
#include <unistd.h> //TEMP - needed for chdir
#endif

using namespace Tahoe;

int main (int argc, char *argv [])
{
#if defined(__MWERKS__) && defined(__MACH__)
	if (chdir("/Volumes/Uster/USERS/tahoe/bin") != 0) cout << " chdir failed" << endl;
	char cwd[255];
	if (getcwd(cwd, 255)) cout << " cwd: " << cwd << endl;
#endif

  MakeCSE_ExecutionT exe;
  sArrayT lineoptions (argc);
  for (int i = 0; i < lineoptions.Length(); i++)
    lineoptions[i] = argv[i];
  exe.Run (lineoptions);
  cout << "\n Program Complete\n\n" << endl;
  return 1;
}
