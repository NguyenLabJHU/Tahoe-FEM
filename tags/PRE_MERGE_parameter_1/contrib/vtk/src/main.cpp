/* $Id: main.cpp,v 1.8 2002-07-02 21:23:00 cjkimme Exp $ */

#include "VTKConsoleT.h"
#include "iConsoleT.h"

using namespace Tahoe;

int main (int argc, char* argv[])
{
  /* list of command-line arguments */
  ArrayT<StringT> arguments(argc);
  for (int i = 0; i < arguments.Length(); i++)
	arguments[i] = argv[i];

  /* construct VTK console object */
  VTKConsoleT vtk_console(arguments);

  /* open console */
  iConsoleT("vtk_console.log", vtk_console, &arguments);


}


