// $Id: main.cpp,v 1.9 2002-11-05 13:26:26 sawimme Exp $
// created: 6 Oct 1999 by S. A. Wimmer
// program reads input file, runs MakeCSE, writes output file

#include "ExceptionT.h"
#include "MakeCSE_ExecutionT.h"
#include "sArrayT.h"

using namespace Tahoe;

int main (int argc, char *argv [])
{
  MakeCSE_ExecutionT exe;
  sArrayT lineoptions (argc);
  for (int i = 0; i < lineoptions.Length(); i++)
    lineoptions[i] = argv[i];
  exe.Run (lineoptions);
  cout << "\n Program Complete\n\n" << endl;
  return 1;
}
